//
//  Created by VT on 08.11.25.
//  © 2025 Volker Thieme. All rights reserved.
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
//  

/// Factory and registry for invoking concrete maximum likelihood fitters.

import SwiftyStatsPrelude

/// Maximum Likelihood Estimation (MLE) namespace for continuous and discrete distributions.
///
/// Use static `fit…` methods on this enum to obtain maximum likelihood estimates
/// for parameters of common probability distributions. When a closed-form (analytic)
/// solution exists, the corresponding method computes it directly. Otherwise, a
/// numerical optimizer is used to maximize the sample log-likelihood.
///
/// - Generic Parameter: `T`
///   The floating-point scalar type used for data and parameters. Must conform to
///   ``RealLike``.
///
/// - Design:
///   - This type is a namespace (an uninstantiable enum) that groups related
///     MLE entry points under a single symbol.
///   - Numerical fits internally construct a solver problem and solve it with a
///     configurable optimizer. Parameter domains and starting values are expressed
///     via internal parameter specifications.
///   - Many methods can optionally estimate an approximate covariance matrix for
///     the parameters at the solution (for example, via the observed information).
///
/// - Return Value:
///   Unless otherwise documented, each `fit…` method returns an ``MLEResult``
///   whose fields include:
///   - `thetaHat`: The MLE parameter vector, ordered as documented by the specific
///     distribution’s `fit…` method.
///   - `logLik`: The log-likelihood evaluated at `thetaHat`.
///   - `iterations`: The number of iterations performed by the optimizer
///     (0 for analytic fits).
///   - `converged`: Whether the optimizer (if used) reported convergence
///     (always `true` for analytic fits).
///   - `nEval`: The number of objective evaluations
///     (0 for analytic fits).
///   - `cov`: An optional covariance estimate in parameter (θ) space, if requested.
///   - `robustCov`: An optional robust (sandwich) covariance when score information is available.
///
/// - Requirements:
///   - All `fit…` methods expect a non-empty `data` array.
///   - Inputs must respect the distribution’s support (for example, strictly positive
///     data for gamma shape/scale fits). Methods commonly validate domains via preconditions.
///   - For numerical fits, parameter constraints are enforced through
///     internal parameter specs and corresponding smooth transforms to an unconstrained space.
///   - The order and parameterization of `thetaHat` are documented per `fit…` method.
///
/// - Discussion:
///   Numerical routines assemble an internal problem from:
///   - the sample data,
///   - a per-observation log-density (and optionally its gradient),
///   - and a set of parameter specifications (constraints, initial values, and u-space steps).
///
///   The problem is then solved by the selected optimizer (for example, Nelder–Mead or L-BFGS).
///   Multi-start strategies, random restarts, and finite-difference settings can be configured
///   through the options type used by your optimizer. When enabled, an approximate covariance
///   matrix is computed at the optimum; a robust covariance may also be provided when score
///   information is available.
///
/// - See Also:
///   - ``RealLike``
///   - ``MLEResult``
///   - ``NelderMeadOptions``
///   - ``OptimizerKind``
///
/// - Note:
///   The concrete `fit…` methods are defined in extensions to this namespace in other files.
///   Each method documents:
///   - the parameterization used,
///   - the order of parameters in `thetaHat`,
///   - any constraints or fixed quantities,
///   - whether the solution is analytic or numerical,
///   - and whether covariance estimates are available.
///
/// - Warning:
///   This enum is not intended to be instantiated.
public enum MLEFitter<T: RealLike> {
}
