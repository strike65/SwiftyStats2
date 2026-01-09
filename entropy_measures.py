#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility script that reproduces the public entropy measures exposed by
`SSExamine` for a fixed numeric series and a German reference text. The
script relies on NumPy only, keeping the implementation self-contained in
the virtual environment described in the project instructions.
"""

from __future__ import annotations

import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


EULER_GAMMA = 0.5772156649015328606

NumericArray = np.ndarray


def probability_mass(values: Sequence) -> Tuple[Counter, np.ndarray, List]:
    """Return Counter, probability vector, and support list for discrete data."""
    counts = Counter(values)
    total = float(sum(counts.values()))
    if total == 0:
        raise ValueError("Cannot compute probabilities of an empty sequence.")
    probs = np.array([count / total for count in counts.values()], dtype=float)
    support = list(counts.keys())
    return counts, probs, support


def log_base(x: np.ndarray, base: float) -> np.ndarray:
    return np.log(x) / np.log(base)


def shannon_entropy(probs: np.ndarray, base: float = 2.0) -> float:
    positive = probs[probs > 0]
    if positive.size == 0:
        return float("nan")
    return float(-np.sum(positive * log_base(positive, base)))


def renyi_entropy(probs: np.ndarray, alpha: float, base: float = 2.0) -> float:
    positive = probs[probs > 0]
    if positive.size == 0:
        return float("nan")
    if alpha == 1.0:
        return shannon_entropy(positive, base=base)
    if alpha == 0.0:
        return float(log_base(len(positive), base))
    if math.isinf(alpha):
        return min_entropy(positive, base=base)
    return float(log_base(np.sum(positive ** alpha), base) / (1.0 - alpha))


def hartley_entropy(probs: np.ndarray, base: float = 2.0) -> float:
    positive = probs[probs > 0]
    return float(log_base(len(positive), base))


def collision_entropy(probs: np.ndarray, base: float = 2.0) -> float:
    positive = probs[probs > 0]
    return float(-log_base(np.sum(positive ** 2), base))


def min_entropy(probs: np.ndarray, base: float = 2.0) -> float:
    positive = probs[probs > 0]
    return float(-log_base(np.max(positive), base))


def tsallis_entropy(probs: np.ndarray, q: float) -> float:
    if q <= 0.0:
        raise ValueError("Tsallis entropy requires q > 0.")
    positive = probs[probs > 0]
    if math.isclose(q, 1.0):
        return float(-np.sum(positive * np.log(positive)))
    return float((1.0 - np.sum(positive ** q)) / (q - 1.0))


def prepare_series(series: Sequence[float]) -> NumericArray:
    data = np.asarray(series, dtype=float)
    if data.ndim != 1:
        raise ValueError("Only 1D series are supported.")
    return data


def approximate_entropy(series: Sequence[float], m: int = 2, r_ratio: float = 0.2) -> float:
    data = prepare_series(series)
    n = data.size
    if n <= m + 1:
        return float("nan")
    std = np.std(data, ddof=0)
    r = r_ratio * std if std > 0 else r_ratio * 1e-12

    def _phi(pattern_length: int) -> float:
        template_count = n - pattern_length + 1
        templates = np.array([data[i : i + pattern_length] for i in range(template_count)])
        distances = np.max(np.abs(templates[:, None, :] - templates[None, :, :]), axis=2)
        comparisons = np.sum(distances <= r, axis=0) / template_count
        safe = np.where(comparisons > 0, comparisons, 1e-12)
        return float(np.mean(np.log(safe)))

    return float(_phi(m) - _phi(m + 1))


def sample_entropy(series: Sequence[float], m: int = 2, r_ratio: float = 0.2) -> float:
    data = prepare_series(series)
    n = data.size
    if n <= m + 1:
        return float("nan")
    std = np.std(data, ddof=0)
    r = r_ratio * std if std > 0 else r_ratio * 1e-12

    def _count(pattern_length: int) -> float:
        template_count = n - pattern_length + 1
        templates = np.array([data[i : i + pattern_length] for i in range(template_count)])
        distances = np.max(np.abs(templates[:, None, :] - templates[None, :, :]), axis=2)
        matches = (distances <= r).astype(float)
        np.fill_diagonal(matches, 0.0)
        return float(np.sum(matches))

    b = _count(m)
    a = _count(m + 1)
    if b == 0.0 or a == 0.0:
        return float("inf")
    return float(-math.log(a / b))


def permutation_entropy(
    series: Sequence[float],
    order: int = 3,
    delay: int = 1,
    base: float = 2.0,
    normalized: bool = True,
) -> float:
    data = prepare_series(series)
    n = data.size
    window = order * delay
    if n < window or order < 2:
        return float("nan")

    patterns: List[Tuple[int, ...]] = []
    for start in range(n - window + 1):
        subseq = data[start : start + window : delay]
        ranks = tuple(np.argsort(np.argsort(subseq, kind="mergesort")))
        patterns.append(ranks)

    counts, probs, _ = probability_mass(patterns)
    h = shannon_entropy(probs, base=base)
    if not normalized:
        return h
    norm = math.log(math.factorial(order), base)
    return float(h / norm if norm > 0 else float("nan"))


def digamma_integer(n: int) -> float:
    if n <= 0:
        raise ValueError("Digamma is undefined for n <= 0.")
    return -EULER_GAMMA + sum(1.0 / k for k in range(1, n))


def kozachenko_leonenko_entropy(series: Sequence[float], k: int = 3, base: float = 2.0) -> float:
    data = prepare_series(series)
    n = data.size
    if n <= k:
        return float("nan")
    distances = np.abs(data[:, None] - data[None, :])
    np.fill_diagonal(distances, np.inf)
    kth = np.partition(distances, k - 1, axis=1)[:, k - 1]
    positive = np.where(kth > 0, kth, np.min(kth[kth > 0]) * 0.5 if np.any(kth > 0) else 1e-12)
    logs = np.log(2.0 * positive)
    h_nat = digamma_integer(n) - digamma_integer(k) + float(np.mean(logs))
    return float(h_nat / np.log(base))


def vasicek_entropy(series: Sequence[float], m: int | None = None, base: float = 2.0) -> float:
    data = np.sort(prepare_series(series))
    n = data.size
    if n < 3:
        return float("nan")
    if m is None:
        m = max(1, int(math.floor(math.sqrt(n))))
    m = max(1, min(m, (n - 1) // 2))
    diffs = data[2 * m :] - data[: n - 2 * m]
    factor = n / (2.0 * m)
    logs = np.log(factor * diffs)
    h_nat = np.sum(logs) / (n - 2.0 * m)
    return float(h_nat / np.log(base))


def conditional_entropy(values: Sequence, condition, base: float = 2.0) -> float:
    subset = [value for value in values if condition(value)]
    if not subset:
        return float("nan")
    _, probs, _ = probability_mass(subset)
    return shannon_entropy(probs, base=base)


def mutual_information(labels_x: Sequence, labels_y: Sequence, base: float = 2.0) -> float:
    if len(labels_x) != len(labels_y):
        raise ValueError("Mutual information requires equally sized label sequences.")
    total = len(labels_x)
    joint_counts = Counter(zip(labels_x, labels_y))
    x_counts = Counter(labels_x)
    y_counts = Counter(labels_y)
    mi = 0.0
    for (x, y), joint in joint_counts.items():
        p_xy = joint / total
        p_x = x_counts[x] / total
        p_y = y_counts[y] / total
        mi += p_xy * math.log(p_xy / (p_x * p_y))
    return float(mi / math.log(base))


def ensure_support_overlap(*counters: Counter) -> List:
    support = set()
    for counter in counters:
        support |= set(counter.keys())
    return list(support)


def cross_entropy(p_counts: Counter, q_counts: Counter, base: float = 2.0) -> float:
    support = ensure_support_overlap(p_counts, q_counts)
    total_p = sum(p_counts.values())
    total_q = sum(q_counts.values())
    ce = 0.0
    for symbol in support:
        p = p_counts.get(symbol, 0.0) / total_p
        q = q_counts.get(symbol, 0.0) / total_q
        if p == 0.0:
            continue
        if q == 0.0:
            return float("inf")
        ce -= p * math.log(q)
    return float(ce / math.log(base))


def kullback_leibler(p_counts: Counter, q_counts: Counter, base: float = 2.0) -> float:
    support = ensure_support_overlap(p_counts, q_counts)
    total_p = sum(p_counts.values())
    total_q = sum(q_counts.values())
    divergence = 0.0
    for symbol in support:
        p = p_counts.get(symbol, 0.0) / total_p
        q = q_counts.get(symbol, 0.0) / total_q
        if p == 0.0:
            continue
        if q == 0.0:
            return float("inf")
        divergence += p * math.log(p / q)
    return float(divergence / math.log(base))


def coarse_grain(series: Sequence[float], scale: int) -> NumericArray:
    data = prepare_series(series)
    n = data.size
    if scale <= 1:
        return data.copy()
    trimmed = n - (n % scale)
    reshaped = data[:trimmed].reshape(-1, scale)
    return np.mean(reshaped, axis=1)


def multiscale_entropy(
    series: Sequence[float],
    scales: Sequence[int] = (1, 2, 3),
    m: int = 2,
    r_ratio: float = 0.2,
) -> float:
    entropies = []
    for scale in scales:
        coarse = coarse_grain(series, scale)
        value = sample_entropy(coarse, m=m, r_ratio=r_ratio)
        if math.isfinite(value):
            entropies.append(value)
    if not entropies:
        return float("nan")
    return float(np.mean(entropies))


def compute_profile(probs: np.ndarray, base: float = 2.0) -> Dict[str, float]:
    profile = {
        "shannon_bits": shannon_entropy(probs, base=base),
        "shannon_nats": shannon_entropy(probs, base=math.e),
        "renyi_alpha_1.5_bits": renyi_entropy(probs, alpha=1.5, base=base),
        "hartley_bits": hartley_entropy(probs, base=base),
        "collision_bits": collision_entropy(probs, base=base),
        "min_bits": min_entropy(probs, base=base),
        "tsallis_q1.5_nats": tsallis_entropy(probs, q=1.5),
    }
    return {key: float(value) for key, value in profile.items()}


def complexity_measures(
    probs: np.ndarray,
    series: Sequence[float] | None = None,
    base: float = 2.0,
) -> Dict[str, float]:
    shannon_nats = shannon_entropy(probs, base=math.e)
    support = np.sum(probs > 0)
    normalized = float(shannon_nats / math.log(support)) if support > 1 else float("nan")
    measures: Dict[str, float] = {
        "normalized_entropy": normalized,
        "statistical_complexity": float(normalized * (1.0 - normalized)) if math.isfinite(normalized) else float("nan"),
    }
    if series is not None:
        mse = multiscale_entropy(series)
        measures["multiscale_entropy"] = mse if math.isfinite(mse) else float("nan")
    else:
        measures["multiscale_entropy"] = float("nan")
    return measures


def sanitize(value):
    if isinstance(value, dict):
        return {k: sanitize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [sanitize(v) for v in value]
    if isinstance(value, (float, np.floating)):
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        if math.isnan(value):
            return "nan"
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    return value


def main() -> None:
    numbers = np.array(
        [
            -0.25,
            0.68,
            0.94,
            1.15,
            1.20,
            1.26,
            1.26,
            1.34,
            1.38,
            1.43,
            1.49,
            1.49,
            1.55,
            1.56,
            1.58,
            1.65,
            1.69,
            1.70,
            1.76,
            1.77,
            1.81,
            1.91,
            1.94,
            1.96,
            1.99,
            2.06,
            2.09,
            2.10,
            2.14,
            2.15,
            2.23,
            2.24,
            2.26,
            2.35,
            2.37,
            2.40,
            2.47,
            2.54,
            2.62,
            2.64,
            2.90,
            2.92,
            2.92,
            2.93,
            3.21,
            3.26,
            3.30,
            3.59,
            3.68,
            4.30,
            4.64,
            5.34,
            5.42,
            6.01,
        ]
    )

    text = (
        "Der alte Ortskern liegt am Schnittpunkt der beiden Hauptstraßen "
        "Kelkheim–Ruppertshain und Königstein–Eppstein. Die erstere bildet im Bereich "
        "des alten Dorfs die Hauptstraße des Ortes, die Langstraße, von der nach beiden "
        "Seiten kurze Gassen abzweigen. Durch einen Graben und eine daran verlaufende "
        "dichte Hecke, die Straße Haingraben erinnert daran, sowie den Lauf des Fischbachs "
        "war der Ort auf einfache Weise vor unliebsamen Besuchern geschützt. Es existierten "
        "auch zwei Wachtürme nebst Toren, welche im Jahre 1348 erstmals erwähnt wurden. "
        "Der untere Turm befand sich ungefähr dort, wo jetzt der Hanseklingerbrunnen steht "
        "(Langstraße, Ecke Kirchgasse). Der obere Turm stand unweit der Langstraße, Ecke "
        "Eppsteiner Straße. Der obere Turm wurde 1818 verkauft und noch im gleichen Jahr abgerissen."
    )

    # Discrete distributions for numbers (exact values) and words (lower-cased tokens).
    number_counts, number_probs, number_support = probability_mass(numbers)
    words = re.findall(r"[\\wäöüÄÖÜß]+", text.lower())
    word_counts, word_probs, word_support = probability_mass(words)

    # Core discrete entropy measures for the numeric series and the word sequence.
    number_measures = {
        "shannon_entropy_bits": shannon_entropy(number_probs),
        "renyi_entropy_alpha_1.5_bits": renyi_entropy(number_probs, alpha=1.5),
        "hartley_entropy_bits": hartley_entropy(number_probs),
        "collision_entropy_bits": collision_entropy(number_probs),
        "min_entropy_bits": min_entropy(number_probs),
        "tsallis_entropy_q1.5_nats": tsallis_entropy(number_probs, q=1.5),
        "approximate_entropy": approximate_entropy(numbers),
        "sample_entropy": sample_entropy(numbers),
        "permutation_entropy_order3_delay1": permutation_entropy(numbers, order=3, delay=1),
        "kozachenko_leonenko_entropy_bits": kozachenko_leonenko_entropy(numbers),
        "vasicek_entropy_bits": vasicek_entropy(numbers),
        "conditional_entropy_high_bits": conditional_entropy(numbers, lambda x: x >= np.median(numbers)),
    }

    # Mutual information from two derived labelings on the numeric series.
    median_value = float(np.median(numbers))
    number_labels_a = ["high" if value >= median_value else "low" for value in numbers]
    number_labels_b = ["even_index" if idx % 2 == 0 else "odd_index" for idx in range(len(numbers))]
    number_measures["mutual_information_bits"] = mutual_information(number_labels_a, number_labels_b)

    number_measures["entropy_profile"] = compute_profile(number_probs)
    number_measures["complexity_measures"] = complexity_measures(number_probs, series=numbers)

    # Discrete entropy measures for the text.
    word_measures = {
        "shannon_entropy_bits": shannon_entropy(word_probs),
        "renyi_entropy_alpha_1.5_bits": renyi_entropy(word_probs, alpha=1.5),
        "hartley_entropy_bits": hartley_entropy(word_probs),
        "collision_entropy_bits": collision_entropy(word_probs),
        "min_entropy_bits": min_entropy(word_probs),
        "tsallis_entropy_q1.5_nats": tsallis_entropy(word_probs, q=1.5),
        "conditional_entropy_long_words_bits": conditional_entropy(words, lambda w: len(w) > np.median([len(word) for word in words])),
    }

    vowel_set = set("aeiouäöü")
    word_labels_a = ["starts_with_vowel" if word[0] in vowel_set else "starts_with_consonant" for word in words]
    word_labels_b = ["short" if len(word) <= 4 else "long" for word in words]
    word_measures["mutual_information_bits"] = mutual_information(word_labels_a, word_labels_b)
    word_measures["entropy_profile"] = compute_profile(word_probs)
    word_measures["complexity_measures"] = complexity_measures(word_probs, series=None)

    # Shared categories for divergence-style reports.
    # Map numeric values and words onto a common binary alphabet to allow KL/Cross computations.
    binary_map = ("A", "B")
    numeric_binary = [binary_map[0] if value < median_value else binary_map[1] for value in numbers]
    median_length = np.median([len(word) for word in words])
    word_binary = [binary_map[0] if len(word) <= median_length else binary_map[1] for word in words]
    numeric_binary_counts, _, _ = probability_mass(numeric_binary)
    word_binary_counts, _, _ = probability_mass(word_binary)

    cross_measures = {
        "cross_entropy_bits": cross_entropy(numeric_binary_counts, word_binary_counts),
        "kullback_leibler_bits": kullback_leibler(numeric_binary_counts, word_binary_counts),
    }

    results = {
        "numbers": sanitize(number_measures),
        "text": sanitize(word_measures),
        "binary_comparison": sanitize(cross_measures),
    }

    output_path = Path("entropy_results.json")
    output_path.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
