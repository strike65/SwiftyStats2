//
//  Created by VT on 30.10.25.
//  Copyright (c) 2025 Volker Thieme. All rights reserved.
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

// This is a testing file - use at your own risk. No structure, no doc
import SwiftyStats
import SwiftyBoost
import Foundation

// Data
let rosnerData: [Double] = [-0.25, 0.68, 0.94, 1.15, 1.20, 1.26, 1.26,
                             1.34, 1.38, 1.43, 1.49, 1.49, 1.55, 1.56,
                             1.58, 1.65, 1.69, 1.70, 1.76, 1.77, 1.81,
                             1.91, 1.94, 1.96, 1.99, 2.06, 2.09, 2.10,
                             2.14, 2.15, 2.23, 2.24, 2.26, 2.35, 2.37,
                             2.40, 2.47, 2.54, 2.62, 2.64, 2.90, 2.92,
                             2.92, 2.93, 3.21, 3.26, 3.30, 3.59, 3.68,
                             4.30, 4.64, 5.34, 5.42, 6.01]
let chars = "Initium mihi operis Servius Galba iterum Titus Vinius consules erunt. nam post conditam urbem octingentos et viginti prioris aevi annos multi auctores rettulerunt, dum res populi Romani memorabantur pari eloquentia ac libertate: postquam bellatum apud Actium atque omnem potentiam ad unum conferri pacis interfuit, magna illa ingenia cessere; simul veritas pluribus modis infracta, primum inscitia rei publicae ut alienae, mox libidine adsentandi aut rursus odio adversus dominantis: ita neutris cura posteritatis inter infensos vel obnoxios. sed ambitionem scriptoris facile averseris, obtrectatio et livor pronis auribus accipiuntur; quippe adulationi foedum crimen servitutis, malignitati falsa species libertatis inest. mihi Galba Otho Vitellius nec beneficio nec iniuria cogniti. dignitatem nostram a Vespasiano inchoatam, a Tito auctam, a Domitiano longius provectam non abnuerim: sed incorruptam fidem professis neque amore quisquam et sine odio dicendus est. quod si vita suppeditet, principatum divi Nervae et imperium Traiani, uberiorem securioremque materiam, senectuti seposui, rara temporum felicitate ubi sentire quae velis et quae sentias dicere licet.[2] Opus adgredior opimum casibus, atrox proeliis, discors seditionibus, ipsa etiam pace saevum. quattuor principes ferro interempti: trina bella civilia, plura externa ac plerumque permixta: prosperae in Oriente, adversae in Occidente res: turbatum Illyricum, Galliae nutantes, perdomita Britannia et statim omissa: coortae in nos Sarmatarum ac Sueborum gentes, nobilitatus cladibus mutuis Dacus, mota prope etiam Parthorum arma falsi Neronis ludibrio. iam vero Italia novis cladibus vel post longam saeculorum seriem repetitis adflicta. haustae aut obrutae urbes, fecundissima Campaniae ora; et urbs incendiis vastata, consumptis, antiquissimis delubris, ipso Capitolio civium manibus incenso. pollutae caerimoniae, magna adulteria: plenum exiliimare, infecti caedibus scopuli. atrocius in urbe saevitum: nobilitas, opes, omissi gestique honores pro crimine et ob virtutes certissimum exitium. nec minus praemia delatorum invisa quam scelera, cum alii sacerdotia et consulatus ut spolia adepti, procurationes alii et interiorem potentiam, agerent verterent cuncta odio et terrore. corrupti in dominos servi, in patronos liberti; et quibus deerat inimicus per amicos oppressi.[3] Non tamen adeo virtutum sterile saeculum ut non et bona exempla prodiderit. comitatae profugos liberos matres, secutae maritos in exilia coniuges: propinqui audentes, constantes generi, contumax etiam adversus tormenta servorum fides; supremae clarorum virorum necessitates fortiter toleratae et laudatis antiquorum mortibus pares exitus. praeter multiplicis rerum humanarum casus caelo terraque prodigia et fulminum monitus et futurorum praesagia, laeta tristia, ambigua manifesta; nec enim umquam atrocioribus populi Romani cladibus magisve iustis indiciis adprobatum est non esse curae deis securitatem nostram, esse ultionem.[4] Ceterum antequam destinata componam, repetendum videtur qualis status urbis, quae mens exercituum, quis habitus provinciarum, quid in toto terrarum orbe validum, quid aegrum fuerit, ut non modo casus eventusque rerum, qui plerumque fortuiti sunt, sed ratio etiam causaeque noscantur. finis Neronis ut laetus primo gaudentium impetu fuerat, ita varios motus animorum non modo in urbe apud patres aut populum aut urbanum militem, sed omnis legiones ducesque conciverat, evulgato imperii arcano posse principem alibi quam Romae fieri. sed patres laeti, usurpata statim libertate licentius ut erga principem novum et absentem; primores equitum proximi gaudio patrum; pars populi integra et magnis domibus adnexa, clientes libertique damnatorum et exulum in spem erecti: plebs sordida et circo ac theatris sueta, simul deterrimi servorum, aut qui adesis bonis per dedecus Neronis alebantur, maesti et rumorum avidi.[5] Miles urbanus longo Caesarum sacramento imbutus et ad destituendum Neronem arte magis et impulsu quam suo ingenio traductus, postquam neque dari donativum sub nomine Galbae promissum neque magnis meritis ac praemiis eundem in pace quem in bello locum praeventamque gratiam intellegit apud principem a legionibus factum, pronus ad novas res scelere insuper Nymphidii Sabini praefecti imperium sibi molientis agitatur. et Nymphidius quidem in ipso conatu oppressus, set quamvis capite defectionis ablato manebat plerisque militum conscientia, nec deerant sermones senium atque avaritiam Galbae increpantium. laudata olim et militari fama celebrata severitas eius angebat aspernantis veterem disciplinam atque ita quattuordecim annis a Nerone adsuefactos ut haud minus vitia principum amarent quam olim virtutes verebantur. accessit Galbae vox pro re publica honesta, ipsi anceps, legi a se militem, non emi; nec enim ad hanc formam cetera erant."


let batch1: [Double] = [608.781, 689.556, 618.134, 680.203, 726.232, 518.655, 740.447, 666.83, 710.272, 751.669, 697.979, 708.583, 624.972, 695.07, 769.391, 720.186, 723.657, 703.7, 697.626, 714.98, 657.712, 609.989, 650.771, 707.977, 712.199, 709.631, 703.16, 744.822, 719.217, 619.137, 753.333, 677.933, 735.919, 695.274, 504.167, 693.333, 625, 596.667, 640.898, 720.506, 700.748, 691.604, 636.738, 731.667, 635.079, 716.926, 759.581, 673.903, 736.648, 675.957, 729.23, 697.239, 728.499, 797.662, 668.53, 815.754, 777.392, 712.14, 663.622, 684.181, 629.012, 640.193, 644.156, 642.469, 639.09, 439.418, 614.664, 537.161, 656.773, 659.534, 695.278, 734.04, 687.665, 710.858, 701.716, 382.133, 719.744, 756.82, 690.978, 670.864, 670.308, 660.062, 790.382, 714.75, 716.959, 603.363, 713.796, 444.963, 723.276, 745.527, 778.333, 723.349, 708.229, 681.667, 566.085, 687.448, 597.5, 637.41, 755.864, 692.945, 766.532, 725.663, 698.818, 760, 775.272, 708.885, 727.201, 642.56, 690.773, 688.333, 743.973, 682.461, 761.43, 691.542, 643.392, 697.075, 708.229, 746.467, 744.819, 655.029, 715.224, 614.417, 761.363, 716.106, 659.502, 730.781, 546.928, 734.203, 682.051, 701.341, 759.729, 689.942, 769.424, 715.286, 776.197, 547.099, 619.942, 696.046, 573.109, 638.794, 708.193, 502.825, 632.633, 683.382, 684.812, 738.161, 671.492, 709.771, 685.199, 624.973, 757.363, 633.417, 658.754, 664.666, 663.009, 773.226, 708.261, 739.086, 667.786, 674.481, 695.688, 588.288, 545.61, 752.305, 684.523, 717.159, 721.343, 750.623, 776.488, 750.623, 600.84, 686.196, 687.87, 725.527, 658.796, 690.38, 737.144, 663.851, 766.63, 625.922, 694.43, 730.217, 700.77, 722.242, 763.828, 695.668, 688.887, 531.021, 698.915, 735.905, 732.039, 751.832, 618.663, 744.845, 690.826, 666.893, 759.86, 683.752, 729.591, 730.706, 763.124, 724.193, 630.352, 750.338, 752.417, 707.899, 715.582, 728.746, 591.193, 592.252, 740.833, 786.367, 712.386, 738.333, 741.48, 729.167, 795.833, 723.502, 718.333, 768.08, 747.5, 775, 760.599, 758.333, 682.5, 658.116, 738.213, 681.236, 704.904, 693.623, 624.993, 700.228, 611.874, 579.167, 720.872, 690.32, 677.933, 674.6, 611.999, 530.68]
let batch2: [Double] = [569.67, 747.541, 612.182, 607.766, 605.38, 589.226, 588.375, 531.384, 633.417, 619.06, 632.447, 624.256, 575.143, 549.278, 624.972, 587.695, 569.207, 613.257, 565.737, 662.131, 543.177, 512.394, 611.19, 659.982, 569.245, 725.792, 608.96, 586.06, 617.441, 592.845, 631.754, 588.113, 555.724, 702.411, 631.754, 698.254, 616.791, 551.953, 636.738, 571.551, 521.667, 587.451, 700.422, 595.819, 534.236, 606.188, 575.303, 590.628, 729.314, 619.313, 624.234, 651.304, 724.175, 583.034, 620.227, 584.861, 565.391, 622.506, 628.336, 587.145, 584.319, 538.239, 538.097, 595.686, 648.935, 583.827, 534.905, 569.858, 617.246, 610.337, 584.192, 598.853, 554.774, 605.694, 627.516, 574.522, 582.682, 563.872, 715.962, 616.43, 778.011, 604.255, 571.906, 625.925, 682.426, 707.604, 617.4, 689.576, 676.678, 563.29, 581.879, 447.701, 557.772, 593.537, 632.585, 671.35, 569.53, 581.667, 643.449, 581.593, 494.122, 620.948, 615.903, 606.667, 579.167, 662.51, 436.237, 644.223, 586.035, 620.833, 652.535, 593.516, 587.451, 570.964, 645.192, 540.079, 707.117, 621.779, 585.777, 703.98, 698.237, 757.12, 621.751, 472.125, 612.7, 583.17, 599.771, 549.227, 605.453, 569.599, 637.233, 621.774, 558.041, 583.17, 345.294, 570.999, 603.232, 595.335, 581.047, 455.878, 627.88, 464.085, 596.129, 640.371, 621.471, 612.727, 606.46, 571.76, 599.304, 579.459, 761.511, 566.969, 654.397, 611.719, 577.409, 576.731, 617.441, 577.409, 548.957, 623.315, 621.761, 553.978, 657.157, 610.882, 552.304, 545.303, 651.934, 635.24, 641.083, 645.321, 566.127, 647.844, 554.815, 620.087, 711.301, 644.355, 713.812, 696.707, 589.453, 634.468, 599.751, 624.542, 723.505, 674.717, 608.539, 612.135, 591.935, 676.656, 647.323, 811.97, 603.883, 608.643, 630.778, 623.063, 472.463, 645.932, 577.176, 567.53, 821.654, 684.49, 600.427, 686.023, 628.109, 605.214, 640.26, 700.767, 665.924, 555.926, 543.299, 511.03, 583.994, 611.048, 623.338, 679.585, 665.004, 655.86, 715.711, 611.999, 577.722, 615.129, 540.316, 711.667, 639.167, 549.491, 684.167, 672.153, 594.534, 627.65, 551.87, 594.534, 602.66, 585.45, 555.724, 574.934, 584.625, 555.724, 611.874, 698.254, 748.13, 689.942]
let gear1: [Double] = [1.006, 0.996, 0.998, 1.000, 0.992, 0.993, 1.002, 0.999, 0.994, 1.000]
let gear2: [Double] = [0.998, 1.006, 1.000, 1.002, 0.997, 0.998, 0.996, 1.000, 1.006, 0.988]
let gear3: [Double] = [0.991, 0.987, 0.997, 0.999, 0.995, 0.994, 1.000, 0.999, 0.996, 0.996]
let gear4: [Double] = [1.005, 1.002, 0.994, 1.000, 0.995, 0.994, 0.998, 0.996, 1.002, 0.996]
let gear5: [Double] = [0.998, 0.998, 0.982, 0.990, 1.002, 0.984, 0.996, 0.993, 0.980, 0.996]
let gear6: [Double] = [1.009, 1.013, 1.009, 0.997, 0.988, 1.002, 0.995, 0.998, 0.981, 0.996]
let gear7: [Double] = [0.990, 1.004, 0.996, 1.001, 0.998, 1.000, 1.018, 1.010, 0.996, 1.002]
let gear8: [Double] = [0.998, 1.000, 1.006, 1.000, 1.002, 0.996, 0.998, 0.996, 1.002, 1.006]
let gear9: [Double] = [1.002, 0.998, 0.996, 0.995, 0.996, 1.004, 1.004, 0.998, 0.999, 0.991]
let gear10: [Double] = [0.991, 0.995, 0.984, 0.994, 0.997, 0.997, 0.991, 0.998, 1.004, 0.997]
let globalSeed: UInt64 = 0x5EED_5EED

// MARK: - Entry Point

/// Entry point wiring the enabled demo routines.
private func runAllDemos() {
    demoMLEAll()
//    demoFileManagement()
//    demoBasicsWithSSExamine()
//    demoZarrStatsAndQuantiles()
    demoOutliersAndEntropy()
//    smokeTests()
}

runAllDemos()


// MARK: - Formatting helpers

@inline(__always)
/// Prints a section banner separating demo phases.
///
/// - Parameter title: Label rendered between `===` markers.
func banner(_ title: String) {
    print("\n=== \(title) ===")
}

@inline(__always)
/// Prints a subsection heading for grouped output.
///
/// - Parameter title: Label rendered between `--` markers.
func section(_ title: String) {
    print("\n-- \(title) --")
}

@inline(__always)
/// Key-value printer used by demo sections.
///
/// - Parameters:
///   - key: Label to print.
///   - value: Any value whose textual description is emitted.
func kv(_ key: String, _ value: Any) {
    print("\(key): \(value)")
}

// MARK: - MLE option helpers

/// Hand tuned option sets that balance robustness and runtime for the demo datasets.
private enum DemoMLEPreset {
    case central      // smooth, 2-3 parameters, well-behaved tails
    case heavyTail    // Cauchy/Holtsmark-style likelihoods with flat regions
    case noncentral   // single-parameter noncentral problems (lambda or delta)
}

/// Builds optimizer presets for the demo MLE runs.
///
/// - Parameters:
///   - preset: Desired preset emphasizing central, heavy-tail, or noncentral problems.
///   - optimizer: Optimizer to use (defaults to L-BFGS).
///   - seed: Optional RNG seed for reproducible multistart layouts.
/// - Returns: Tuned `NelderMeadOptions` shared across demo fitters.
private func demoOptions(_ preset: DemoMLEPreset,
                         optimizer: OptimizerKind = .lbfgs,
                         seed: UInt64? = nil) -> MLEOptimizationOpts<Double> {
    var opt = MLEOptimizationOpts<Double>()
    opt.optimizer = optimizer
    opt.computeCovariance = true
    opt.rngSeed = seed
    opt.multiStartDesign = .sobol
    opt.diagnosticsEnabled = true
    switch preset {
    case .central:
        opt.multiStartCount = 6
        opt.randomRestartCount = 2
        opt.initialStepStrategy = .relativeTheta(0.3)
        opt.gradStep = 5e-6
        opt.hessianStep = 1e-3
        opt.relTolLogLik = 1e-7
    case .heavyTail:
        opt.multiStartCount = 8
        opt.randomRestartCount = 4
        opt.initialStepStrategy = .relativeTheta(0.4)
        opt.gradStep = 8e-6
        opt.hessianStep = 2e-3
        opt.boundaryPenaltyWeight = 1e-3
        opt.relTolLogLik = 5e-8
    case .noncentral:
        opt.multiStartCount = 10
        opt.randomRestartCount = 4
        opt.initialStepStrategy = .relativeTheta(0.5)
        opt.gradStep = 2e-5
        opt.hessianStep = 7.5e-4
        opt.relTolLogLik = 1e-8
    }
    return opt
}

// MARK: - Demos

/// Demonstrates SSExamine ingestion, ordering, and frequency rendering for integer and string samples.
private func demoBasicsWithSSExamine() {
    banner("Basics with SSExamine")

    var examine1: SSExamine<Int, Double>
    var examine2: SSExamine<String, Double>

    do {
        section("Ordinal Int sample")
        examine1 = try SSExamine(using: [1,2,3,4,5,43,23,2,3,4,54], levelOfMeasurement: .ordinal, name: "Examine1", characterSet: nil)
        var array = examine1.itemsAsArray(sorted: .original)
        kv("original", array)
        array = examine1.itemsAsArray(sorted: .ascending)
        kv("ascending", array)
        array = examine1.itemsAsArray(sorted: .descending)
        kv("descending", array)
        array = examine1.itemsAsArray(sorted: .none)
        kv("none", array)
        kv("sampleSize", examine1.sampleSize)
        examine1.removeAll()

        section("Nominal String sample (letters only)")
        examine2 = try SSExamine(using: "EeHalllooooem,dmnmdkcnn", levelOfMeasurement: .nominal, name: nil, characterSet: .alphanumerics)
        var array1 = examine2.itemsAsArray(sorted: .original)
        kv("original", array1)
        array1 = examine2.itemsAsArray(sorted: .ascending)
        kv("ascending", array1)
        array1 = examine2.itemsAsArray(sorted: .descending)
        kv("descending", array1)
        array1 = examine2.itemsAsArray(sorted: .none)
        kv("none", array1)
        examine2.remove("e", allOccurences: false)
        if let s = examine2.itemsAsString(delimiter: "", asRow: true, encloseElementsBy: "\"") {
            kv("itemsAsString", s)
        } else {
            kv("itemsAsString", String(describing: examine2.itemsAsString))
        }
    } catch {
        kv("Error", error.localizedDescription)
    }
}

/// Shows descriptive statistics, confidence intervals, and Hyndman-Fan quantiles on a Zarr-derived dataset.
private func demoZarrStatsAndQuantiles() {
    banner("Zarr data: Descriptives and Quantiles")
    let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/Zarr.csv", directoryHint: .notDirectory)
    if let zarrExamine: SSExamine<Double, Double> = .loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
        section("Descriptives")
        kv("mean", zarrExamine.arithmeticMean!)
        kv("sample SD", zarrExamine.sampleStandardDeviation!)
        kv("sample variance", zarrExamine.sampleVariance!)
        kv("mode", zarrExamine.mode!)
        
        section("Student-t CI (alpha = 0.05)")
        let ci: ConfidenceInterval<Double> = try! zarrExamine.studentTCI(alpha: 0.05)!
        print(ci.description)
        
        section("Normal quantile(0.75) reciprocal")
        kv("1 / Phi^{-1}(0.75)", try! SwiftyBoost.Distribution.Normal().quantile(0.75).reciprocal ?? 0)
        
        section("File identifier")
        kv("#fileID", #fileID)
        
        section("Quantiles HF2 default")
        if let q25 = try? zarrExamine.quantile(q: 0.25),
           let q50 = try? zarrExamine.quantile(q: 0.50),
           let q75 = try? zarrExamine.quantile(q: 0.75) {
            kv("Q1", q25)
            kv("Median", q50)
            kv("Q3", q75)
        } else {
            kv("Quantiles", "Unable to compute (HF2).")
        }
        
        section("Quartile convenience (HF2)")
        if let qs = try? zarrExamine.quartile() {
            print(qs.description)
        } else {
            kv("Quartiles", "Unavailable")
        }
        
        section("Quantiles for Hyndman-Fan Types 1...9")
        let hfTypes: [HyndmanFanQuantileType] = [
            .type1, .type2, .type3, .type4, .type5, .type6, .type7, .type8, .type9
        ]
        for t in hfTypes {
            let q1 = try? zarrExamine.quantile(q: 0.25, quantileType: t)
            let q2 = try? zarrExamine.quantile(q: 0.50, quantileType: t)
            let q3 = try? zarrExamine.quantile(q: 0.75, quantileType: t)
            print("[\(t.description)] Q1=\(String(describing: q1)), Q2/Median=\(String(describing: q2)), Q3=\(String(describing: q3))")
            if let qts = try? zarrExamine.quartile(quantileType: t) {
                print("  Quartile summary - q25: \(qts.lower), q50: \(qts.median), q75: \(qts.upper)")
            } else {
                print("  Quartile summary unavailable for \(t.description)")
            }
        }
    }
}

/// Exercises Grubbs and Rosner ESD outlier tests plus entropy estimators for numeric and string data.
private func demoOutliersAndEntropy() {
    banner("Outliers, ESD tests, and Entropy")

    section("Grubbs example data")
    let grubbsData: [Double] = [199.31,199.53,200.19,200.82,201.92,201.95,202.18,245.57]
    let grubbsEx = try! SSExamine<Double,Double>(using: grubbsData, levelOfMeasurement: .ratio, name: "Grubbs test example data", characterSet: nil)
    let g1 = try! Outliers<Double>.grubbsTest(values: grubbsData, mean: grubbsEx.arithmeticMean!, standardDeviation: grubbsEx.sampleStandardDeviation!, alpha: 0.05)!
    kv("Grubbs test", g1)
    var esdTest = grubbsEx.rosnerESDTest(maxOutliers: 5, tails: .bothTails)!
    kv("Rosner ESD", esdTest)

    section("Rosner example data")
    let rosnerEx = try! SSExamine<Double,Double>(using: rosnerData, levelOfMeasurement: .ratio, name: "Rosner test example data", characterSet: nil)
    let g2 = try! Outliers<Double>.grubbsTest(values: rosnerData, mean: rosnerEx.arithmeticMean!, standardDeviation: rosnerEx.sampleStandardDeviation!, alpha: 0.05)!
    kv("Grubbs test", g2)
    kv("Outliers (both tails)", rosnerEx.outliers(alpha: 0.05, maxOutliers: 10, tails: .bothTails)!)
    kv("Outliers (grubbs data, both tails)", grubbsEx.outliers(alpha: 0.05, maxOutliers: 3, tails: .bothTails)!)
    esdTest = rosnerEx.rosnerESDTest(maxOutliers: 10, tails: .bothTails)!
    kv("Rosner ESD", esdTest)

    section("Entropy")
    kv("Shannon (histogram)", rosnerEx.shannonEntropy(type: .histogram, base: Double.exp(1), bins: .auto))
    kv("Shannon (discrete)", rosnerEx.shannonEntropy(type: .discrete, base: Double.exp(1)))
    kv("Shannon (frequency)", rosnerEx.shannonEntropy(type: .frequency, base: Double.exp(1)))
    kv("Shannon (kde,bw=2.0)", rosnerEx.shannonEntropy(type: .kde(bandwidth: 2.0), base: Double.exp(1)))
    kv("Shannon (Miller-Madow, base 2, Sturges)", rosnerEx.shannonEntropy(type: .millerMadow, base: 2, bins: .sturges))
    kv("Renyi alpha=0.05, base 2", rosnerEx.renyiEntropy(alpha: 0.05, base: 2))

    section("String entropy and divergence")
    let string: String = "ABCEDEFS#X"
    let string1 = "ABCEDEFS#X4e"
    let stringEx : SSExamine<String, Double> = try! SSExamine<String, Double>(using: string, levelOfMeasurement: .nominal, name: nil, characterSet: .letters)
    kv("Shannon (histogram)", stringEx.shannonEntropy(type: .histogram, base: Double.exp(1)))
    kv("Shannon (discrete)", stringEx.shannonEntropy(type: .discrete, base: Double.exp(1)))
    kv("Complexity measures", stringEx.complexityMeasures())
    kv("Entropy profile", stringEx.computeEntropyProfile())
    let stringEx1 : SSExamine<String, Double> = try! SSExamine<String, Double>(using: string1, levelOfMeasurement: .nominal, name: nil, characterSet: .letters)
    kv("Kullback-Leibler D(string || string1)", stringEx.kullbackLeiblerDivergence(to: stringEx1))
}

/// Runs a handful of representative MLE fits to illustrate basic solver usage and reporting.
private func demoMLEQuickStart() {
    banner("MLE quick start")
    var result: MLEResult<Double> = MLEFitter<Double>.fitSkewNormal(
        rosnerData,
        optimizer: .lbfgs,
        options: demoOptions(.heavyTail, seed: 0x5EED_5EED)
    )
    print("thetaHat =", result.thetaHat, " logLik =", result.logLik, " iter =", result.iterations, " converged =", result.converged)
    if let cov = result.cov { print("Cov(thetaHat) approx", cov) }

    // Beta sample via RNGSampler.randoms (n = 300)
    let betaData: [Double] = RNGSampler<Double>.randoms(n: 10_000, dist: .beta(a: 22, b: 0.5))
    print(betaData)
    result = MLEFitter<Double>.fitBeta(
        betaData,
        optimizer: .lbfgs,
        options: demoOptions(.central, seed: 0xBADA_55)
    )
    print("thetaHat =", result.thetaHat, " (true [alpha=22.0, beta=0.5])", " logLik =", result.logLik, " iter =", result.iterations, " converged =", result.converged)

    // Gumbel sample via RNGSampler.randoms (n = 300)
    let zr: [Double] = RNGSampler<Double>.randoms(n: 10_000, dist: .gumbel(mu: 3, beta: 5))
    result = MLEFitter<Double>.fitExtremeValue(
        zr,
        optimizer: .lbfgs,
        options: demoOptions(.central, seed: 0xE17E_0003)
    )
    print("#############################\nDistribution: Gumbel\nthetaHat =", result.thetaHat, " (true [mu=3.0, beta=5.0])", " logLik =", result.logLik, " iter =", result.iterations, " converged =", result.converged)
}

/// Executes every available MLE fitter with reproducible multistart options to showcase breadth.
private func demoMLEAll() {
    banner("MLE demos (comprehensive)")
    // Exponential (analytic)
//    do {
//        let rate = 2.5
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-exponential.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let dd = dataEx.itemsAsNumericArray {
//                data = dd
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 10_000, dist: .exponential(rate: rate), seed: globalSeed)
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 10_000, dist: .exponential(rate: rate), seed: globalSeed)
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "exponential")
//        }
//        let r = MLEFitter<Double>.fitExponential(data: data)
//        print("#############################\nDistribution: Exponential(rate)")
//        print("thetaHat =", r.thetaHat, " (true [rate=\(rate)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Rayleigh (sigma): derive from U(0,1) using R = sigma sqrt(-2 ln U)
//    do {
//        let sigma = 1.2
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-rayleigh.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                let u = RNGSampler<Double>.randoms(n: 10_000, dist: .uniform01, seed: globalSeed)
//                data = u.map { sigma * sqrt(-2.0 * log(max($0, .leastNonzeroMagnitude))) }
//            }
//        }
//        else {
//            let u = RNGSampler<Double>.randoms(n: 10_000, dist: .uniform01, seed: globalSeed)
//            data = u.map { sigma * sqrt(-2.0 * log(max($0, .leastNonzeroMagnitude))) }
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "rayleigh")
//        }
//        let r = MLEFitter<Double>.fitRayleigh(data: data)
//        print("#############################\nDistribution: Rayleigh")
//        print("thetaHat =", r.thetaHat, " (true [sigma=\(sigma), sigma^2=\(sigma * sigma)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//
//    }
//
//    // Laplace(mu, b): derive from U(0,1)
//    do {
//        let mu = -10.0, b = 3.0
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-Laplace.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                let u = RNGSampler<Double>.randoms(n: 10_000, dist: .uniform01, seed: globalSeed)
//                let centered = u.map { $0 - 0.5 }
//                data = centered.map { t -> Double in
//                    mu - b * (t >= 0 ? 1 : -1) * log(1 - 2 * abs(t))
//                }
//            }
//        }
//        else {
//            let u = RNGSampler<Double>.randoms(n: 10_000, dist: .uniform01, seed: globalSeed)
//            let centered = u.map { $0 - 0.5 }
//            data = centered.map { t -> Double in
//                mu - b * (t >= 0 ? 1 : -1) * log(1 - 2 * abs(t))
//            }
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "Laplace")
//        }
//        let r = MLEFitter<Double>.fitLaplace(data: data)
//        print("#############################\nDistribution: Laplace")
//        print("thetaHat =", r.thetaHat, " (true [mu=\(mu), b=\(b)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Uniform(a, b): derive from U(0,1)
//    do {
//        let a = -2.0, b = 5.0
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-Uniform.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                let u = RNGSampler<Double>.randoms(n: 10_000, dist: .uniform01, seed: globalSeed)
//                data = u.map { a + (b - a) * $0 }
//            }
//        }
//        else {
//            let u = RNGSampler<Double>.randoms(n: 10_000, dist: .uniform01, seed: globalSeed)
//            data = u.map { a + (b - a) * $0 }
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "Uniform")
//        }
//        let r = MLEFitter<Double>.fitUniform(data: data)
//        print("#############################\nDistribution: Uniform")
//        print("thetaHat =", r.thetaHat, " (true [\(a), \(b)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Wald (Inverse Gaussian): use normals and uniforms via randoms
//    do {
//        let mu = 1.0, lambda = 2.5
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-Wald.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                let z = RNGSampler<Double>.randoms(n: 10_000, dist: .gaussian(mean: 0, standardDeviation: 1), seed: globalSeed)
//                let u = RNGSampler<Double>.randoms(n: 10_000, dist: .uniform01, seed: globalSeed)
//                data = zip(z, u).map { (zVal, uVal) -> Double in
//                    let y = zVal * zVal
//                    let muOverLam = mu / lambda
//                    let term = mu + (mu * mu * y) * 0.5 * (1.0 / lambda) - (muOverLam * 0.5) * sqrt(4 * mu * lambda * y + mu * mu * y * y)
//                    return (uVal <= mu / (mu + term)) ? term : (mu * mu / term)
//                }
//            }
//        }
//        else {
//            let z = RNGSampler<Double>.randoms(n: 10_000, dist: .gaussian(mean: 0, standardDeviation: 1), seed: globalSeed)
//            let u = RNGSampler<Double>.randoms(n: 10_000, dist: .uniform01, seed: globalSeed)
//            data = zip(z, u).map { (zVal, uVal) -> Double in
//                let y = zVal * zVal
//                let muOverLam = mu / lambda
//                let term = mu + (mu * mu * y) * 0.5 * (1.0 / lambda) - (muOverLam * 0.5) * sqrt(4 * mu * lambda * y + mu * mu * y * y)
//                return (uVal <= mu / (mu + term)) ? term : (mu * mu / term)
//            }
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "Wald")
//        }
//        let r = MLEFitter<Double>.fitWald(data: data)
//        print("#############################\nDistribution: Wald (Inverse Gaussian)")
//        print("thetaHat =", r.thetaHat, " (true [mu=\(mu), lambda=\(lambda)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Arcsine(a, b)
//    do {
//        let a = -1.5, b = 2.0
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-ArcSine.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data  = RNGSampler<Double>.randoms(n: 10_000, dist: .arcsine(a: a, b: b))
//            }
//        }
//        else {
//            data  = RNGSampler<Double>.randoms(n: 10_000, dist: .arcsine(a: a, b: b))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "ArcSine")
//        }
//        let r = MLEFitter<Double>.fitArcsine(data: data)
//        print("#############################\nDistribution: Arcsine")
//        print("thetaHat =", r.thetaHat, " (true [a=\(a), b=\(b)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Beta(alpha, beta)
//    do {
//        let a = 2.2, b = 0.6
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-beta.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 10_000, dist: .beta(a: a, b: b))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 10_000, dist: .beta(a: a, b: b))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "beta")
//        }
//        let r = MLEFitter<Double>.fitBeta(
//            data,
//            optimizer: .lbfgs,
//            options: demoOptions(.central, seed: 0x00BA_DA55)
//        )
//        print("#############################\nDistribution: Beta")
//        print("thetaHat =", r.thetaHat, " (true [alpha=\(a), beta=\(b)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Gamma(k, theta)
//    do {
//        let k = 2.0, th = 1.35
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-gamma.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 1000, dist: .gamma(shape: k, scale: th))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 1000, dist: .gamma(shape: k, scale: th))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "gamma")
//        }
//        let r = MLEFitter<Double>.fitGamma(
//            data,
//            optimizer: .lbfgs,
//            options: demoOptions(.central, seed: 0xFACE_B00C)
//        )
//        print("#############################\nDistribution: Gamma")
//        print("thetaHat =", r.thetaHat, " (true [k=\(k), theta=\(th)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Weibull(k, lambda)
//    do {
//        let k = 1.4, lambda = 0.9
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-weibull.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 10_000, dist: .weibull(k: k, lambda: lambda))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 10_000, dist: .weibull(k: k, lambda: lambda))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "weibull")
//        }
//        let r = MLEFitter<Double>.fitWeibull(
//            data,
//            optimizer: .lbfgs,
//            options: demoOptions(.central, seed: 0xDEAD_C0DE)
//        )
//        print("#############################\nDistribution: Weibull")
//        print("thetaHat =", r.thetaHat, " (true [k=\(k), lambda=\(lambda)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Cauchy(mu, gamma)
//    do {
//        let mu = 0.3, gamma = 1.1
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-cauchy.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 10_000, dist: .cauchy(mu: mu, gamma: gamma))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 10_000, dist: .cauchy(mu: mu, gamma: gamma))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "cauchy")
//        }
//        let r = MLEFitter<Double>.fitCauchy(
//            data,
//            optimizer: .lbfgs,
//            options: demoOptions(.heavyTail, seed: 0xC0C0_A)
//        )
//        print("#############################\nDistribution: Cauchy")
//        print("thetaHat =", r.thetaHat, " (true [mu=\(mu), gamma=\(gamma)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Logistic(mu, s)
//    do {
//        let mu = -0.75, s = 0.65
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-logistic.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 10_000, dist: .logistic(mu: mu, s: s))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 10_000, dist: .logistic(mu: mu, s: s))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "logistic")
//        }
//        let r = MLEFitter<Double>.fitLogistic(
//            data,
//            optimizer: .lbfgs,
//            options: demoOptions(.central, seed: 0xC0FF_EE)
//        )
//        print("#############################\nDistribution: Logistic")
//        print("thetaHat =", r.thetaHat, " (true [mu=\(mu), s=\(s)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Student's t (nu only in fitter; generate standard t)
//    do {
//        let nu = 7.5
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-student.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 10_000, dist: .studentT(mu: 0, sigma: 1, nu: nu))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 10_000, dist: .studentT(mu: 0, sigma: 1, nu: nu))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "student")
//        }
//        let r = MLEFitter<Double>.fitStudentsT(
//            data,
//            optimizer: .lbfgs,
//            options: demoOptions(.central, seed: 0x0707_E7E7)
//        )
//        print("#############################\nDistribution: Student's t (nu only)")
//        print("thetaHat =", r.thetaHat, " (true [nu=\(nu)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Fisher F(d1, d2)
//    do {
//        let d1 = 10.0, d2 = 14.0
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-fisher.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 10_000, dist: .f(d1: d1, d2: d2))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 10_000, dist: .f(d1: d1, d2: d2))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "fisher")
//        }
//        let r = MLEFitter<Double>.fitFisherF(
//            data,
//            optimizer: .lbfgs,
//            options: demoOptions(.central, seed: 0x0F17_F17)
//        )
//        print("#############################\nDistribution: Fisher F")
//        print("thetaHat =", r.thetaHat, " (true [d1=\(d1), d2=\(d2)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Inverse-Gamma(alpha, beta): derive from Gamma(alpha,1)
//    do {
//        let alpha = 3.0, beta = 2.0
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-inverseGamma.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                let y = RNGSampler<Double>.randoms(n: 10_000, dist: .gamma(shape: alpha, scale: 1))
//                data = y.map { beta / $0 }
//            }
//        }
//        else {
//            let y = RNGSampler<Double>.randoms(n: 10_000, dist: .gamma(shape: alpha, scale: 1))
//            data = y.map { beta / $0 }
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "inverseGamma")
//        }
//        var opts = demoOptions(.central, seed: 0x01AF_E0)
//        let mean = data.reduce(0.0, +) / Double(data.count)
//        let variance = max(data.reduce(0.0) { $0 + ($1 - mean) * ($1 - mean) } / Double(data.count - 1), 1e-6)
//        let alphaGuess = max(2.05, 2.0 + (mean * mean) / variance)
//        let betaGuess = max(1e-3, mean * (alphaGuess - 1.0))
//        opts.warmStartTheta = [alphaGuess, betaGuess]
//        let r = MLEFitter<Double>.fitInverseGamma(data, optimizer: .lbfgs, options: opts)
//        print("#############################\nDistribution: Inverse-Gamma")
//        print("thetaHat =", r.thetaHat, " (true [alpha=\(alpha), beta=\(beta)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Skew-Normal(xi, omega, alpha)
//    do {
//        let xi = -0.2, omega = 1.1, alpha = 3.0
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-skewNormal.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 10_000, dist: .skewNormal(xi: xi, omega: omega, alpha: alpha))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 10_000, dist: .skewNormal(xi: xi, omega: omega, alpha: alpha))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "skewNormal")
//        }
//        let r = MLEFitter<Double>.fitSkewNormal(
//            data,
//            optimizer: .lbfgs,
//            options: demoOptions(.heavyTail, seed: 0x05AE_E5)
//        )
//        print("#############################\nDistribution: Skew-Normal")
//        print("thetaHat =", r.thetaHat, " (true [xi=\(xi), omega=\(omega), alpha=\(alpha)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Noncentral Chi-Squared(k, lambda)
//    do {
//        let k = 6.0, lambda = 5.0
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-nc-chiSquare.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 10_000, dist: .noncentralChiSquare(k: k, lambda: lambda))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 10_000, dist: .noncentralChiSquare(k: k, lambda: lambda))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "nc-chiSquare")
//        }
//        let r = MLEFitter<Double>.fitNCChiSquared(
//            data,
//            degreesOfFreedom: k,
//            optimizer: .lbfgs,
//            options: demoOptions(.noncentral, seed: 0x0AC2_C2)
//        )
//        print("#############################\nDistribution: Noncentral Chi-Squared")
//        print("thetaHat =", r.thetaHat, " (true [lambda=\(lambda)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Noncentral Student's t(nu, delta) - generate standard (mu=0, sigma=1)
//    do {
//        let nu = 9.0, delta = 5.0
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-nc-student.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 10_000, dist: .noncentralT(mu: 0, sigma: 1, nu: nu, delta: delta))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 10_000, dist: .noncentralT(mu: 0, sigma: 1, nu: nu, delta: delta))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "nc-student")
//        }
//        let r = MLEFitter<Double>.fitNCStudentsT(
//            data,
//            degreesOfFreedom: nu,
//            optimizer: .lbfgs,
//            options: demoOptions(.noncentral, seed: 0x00AC_E7)
//        )
//        print("#############################\nDistribution: Noncentral Student's t")
//        print("thetaHat =", r.thetaHat, " (true [delta=\(delta)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
//
//    // Noncentral Fisher F(d1, d2, lambda): compose from nc-chi^2 and chi^2
//    do {
//        let d1 = 8.0, d2 = 16.0, lambda = 3.5
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-nc-fisher.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                let x1 = RNGSampler<Double>.randoms(n: 10_000, dist: .noncentralChiSquare(k: d1, lambda: lambda))
//                let x2 = RNGSampler<Double>.randoms(n: 10_000, dist: .chiSquare(df: d2))
//                data = zip(x1, x2).map { (a, b) in (a / d1) / (b / d2) }
//            }
//        }
//        else {
//            let x1 = RNGSampler<Double>.randoms(n: 10_000, dist: .noncentralChiSquare(k: d1, lambda: lambda))
//            let x2 = RNGSampler<Double>.randoms(n: 10_000, dist: .chiSquare(df: d2))
//            data = zip(x1, x2).map { (a, b) in (a / d1) / (b / d2) }
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "nc-fisher")
//        }
//        var opt = demoOptions(.noncentral, seed: 0x0ACF_001 ^ 0x0ACF_002)
//        let meanF = max((try? SSExamine<Double, Double>(using: data, levelOfMeasurement: .ratio, name: nil, characterSet: nil).arithmeticMean) ?? 1.0, 1.000001)
//        let lamWarm = max(((d1 * (d2 - 2.0)) / d2) * meanF - d1, 1e-6)
//        opt.warmStartTheta = [lamWarm]
//        
//        // Now fit
//        let r = MLEFitter<Double>.fitNCFisherF(data,df1: d1, df2: d2, optimizer: .lbfgs, options: opt)
//        print("ncF thetaHat =", r.thetaHat, "logLik =", r.logLik, "iters =", r.iterations, "converged =", r.converged)
//        if let sols = r.allSolutions {
//            print("All terminal solutions (thetaHat, logLik):")
//            for (th, ll) in sols { print(th, ll) }
//            print("Unique u-solutions:", r.uniqueSolutionCount ?? 0)
//            if let uniqTheta = r.uniqueSolutionCountTheta {
//                print("Unique theta-solutions:", uniqTheta)
//            }
//            if let cond = r.conditionNumberEstimateHu {
//                let source = r.conditionSource ?? "H_u"
//                print("Condition estimate [\(source)]:", cond)
//            }
//            if let reason = r.convergenceReason {
//                print("Convergence reason:", reason.rawValue)
//            }
//            print("#############################\nDistribution: Noncentral Fisher F")
//            print("thetaHat =", r.thetaHat, " (true [d1=\(d1), d2=\(d2), lambda=\(lambda)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//        }
//    }
//
//    // Noncentral Beta(alpha, beta, lambda): compose from (nc-chi^2, chi^2)
//    do {
//        // Improve estimation stability via: larger n, method-of-moments warm start,
//        // parameter scaling, light bounds, and a few multi-starts.
//        let a = 3.0, b = 2.5, lambda = 1.2
//        let n = 10_000 // larger sample improves identifiability for lambda
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-nc-beta.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true),
//           let d = dataEx.itemsAsNumericArray {
//            data = d
//        } else {
//            let X = RNGSampler<Double>.randoms(n: n, dist: .noncentralChiSquare(k: 2 * a, lambda: 2 * lambda), seed: globalSeed)
//            let Y = RNGSampler<Double>.randoms(n: n, dist: .chiSquare(df: 2 * b), seed: globalSeed ^ 0xA5A5)
//            data = zip(X, Y).map { $0 / ($0 + $1) }
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "nc-beta")
//        }
//
//        // Moment-based warm start from the observed sample
//        let ex: SSExamine<Double, Double> = .init(usingArray: data, name: "nc-beta", characterSet: nil)
//        let m = max(min(ex.arithmeticMean ?? 0.5, 0.99), 0.01)
//        let v = max(min(ex.sampleVariance ?? 0.02, 0.25), 1e-6)
//        // For central Beta, MoM: alpha0 = m*((m*(1-m))/v - 1), beta0 = (1-m)*((m*(1-m))/v - 1)
//        // Use these as a conservative anchor even for noncentral; let optimizer adjust lambda.
//        let t = max((m * (1 - m)) / v - 1.0, 2.05) // ensure > 1 for stability
//        let alpha0 = max(m * t, 1.05)
//        let beta0  = max((1 - m) * t, 1.05)
//        let lambda0 = max(0.2, min(5.0, 2.0 * abs((m - a/(a+b)))) ) // small positive start
//
//        var opts = demoOptions(.noncentral, optimizer: .lbfgs, seed: 0x0ACB_001)
//        // Strengthen search a bit
//        opts.multiStartCount = max(opts.multiStartCount, 8)
//        opts.randomRestartCount = max(opts.randomRestartCount, 3)
//        opts.boundaryPenaltyWeight = max(opts.boundaryPenaltyWeight, 1e-2)
//        // Provide warm start and light parameter scaling
//        opts.warmStartTheta = [alpha0, beta0, lambda0]
//
//        let r = MLEFitter<Double>.fitNCBeta(
//            data,
//            optimizer: .lbfgs,
//            options: opts
//        )
//        print("#############################\nDistribution: Noncentral Beta")
//        print("thetaHat =", r.thetaHat, " (true [alpha=\(a), beta=\(b), lambda=\(lambda)])", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//        if let cov = r.cov { print("Cov(thetaHat) approx", cov) }
//        if let reason = r.convergenceReason { print("Convergence reason:", reason.rawValue) }
//    }
//
//    // Holtsmark
//    do {
//        let mu = -1.5, c = 1.4
//        let data: [Double]
//        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-holtsmark.csv", directoryHint: .notDirectory)
//        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
//            if let d = dataEx.itemsAsNumericArray {
//                data = d
//            }
//            else {
//                data = RNGSampler<Double>.randoms(n: 1000, dist: .holtsmark(mu: mu, c: c ))
//            }
//        }
//        else {
//            data = RNGSampler<Double>.randoms(n: 1000, dist: .holtsmark(mu: mu, c: c ))
//            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "holtsmark")
//        }
//        let r = MLEFitter<Double>.fitHoltsmark(
//            data,
//            optimizer: .lbfgs,
//            options: demoOptions(.heavyTail, seed: 0x0A17_0005)
//        )
//        print("#############################\nDistribution: Holtsmark")
//        print("thetaHat =", r.thetaHat, " (true [mu=\(mu), c=\(c))", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
//    }
    do {
        let a: Double = 0, b: Double = 1
        let data: [Double]
        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/RND-Landau.csv", directoryHint: .notDirectory)
        if let dataEx = SSExamine<Double, Double>.loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
            if let d = dataEx.itemsAsNumericArray {
                data = d
            }
            else {
                data = RNGSampler<Double>.randoms(n: 10_000, dist: .landau(location: a, scale: b))
            }
        }
        else {
            data = RNGSampler<Double>.randoms(n: 10_000, dist: .holtsmark(mu: a, c: b ))
            writeTestValues(examine: SSExamine<Double,Double>.init(usingArray: data, name: nil, characterSet: nil), name: "Landau")
        }
        let r = MLEFitter<Double>.fitLandau(
            data ,
            optimizer: .lbfgs,
            options: demoOptions(.heavyTail, seed: 0x0A17_0005)
        )
        print("#############################\nDistribution: Landau")
        print("thetaHat =", r.thetaHat, " (true [a=\(a), b=\(b))", " logLik =", r.logLik, " iter =", r.iterations, " converged =", r.converged)
    }


}

/// Demonstrates file import/export round-trips for SSExamine in text, JSON, and ZIP formats.
private func demoFileManagement() {
    print("File management demo")
    do {
        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/Zarr.csv", directoryHint: .notDirectory)
        if let examineDoubles: SSExamine<Double, Double> = .loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
            let examineString:SSExamine<String, Double> = try SSExamine<String, Double>.init(using: chars, levelOfMeasurement: .nominal, name: "Tacitus", characterSet: .alphanumerics)
            let sourceDir = URL(filePath: #filePath).deletingLastPathComponent()
            var outURL = sourceDir.appending(path: "Zarr.txt", directoryHint: .notDirectory)
            // Doubles
            let _ = try examineDoubles.saveTo(fileName: outURL.path, atomically: true, overwrite: true, separator: ",", encloseElementsBy: nil, asRow: true, stringEncoding: String.Encoding.utf8, compressAsZip: true, zipEntryFileName: nil)
            outURL = sourceDir.appending(path: "Zarr.txt.zip")
            if let examine2: SSExamine<Double, Double> = try .examine(fromFile: outURL.path, separator: ",", { token in
                guard let t = token?.trimmingCharacters(in: .whitespacesAndNewlines), !t.isEmpty else { return nil }
                return Double(t)
            } )
            {
                print(examineDoubles.isEqual(examine2))
            }
            else {
                print("error")
            }
            // cleanup Zarr.txt.zip
            do { try FileManager.default.removeItem(at: outURL) } catch { print("Cleanup failed (\(outURL.lastPathComponent)): \(error.localizedDescription)") }
            
            outURL = sourceDir.appending(path: "Zarr.json")
            let _ = try examineDoubles.exportJSON(fileName: outURL.path, atomically: true, overwrite: true)
            if let jsonExamine: SSExamine<Double, Double> = try .examine(fromJSON: outURL.path) {
                print(examineDoubles.isEqual(jsonExamine))
            }
            // cleanup Zarr.json
            do { try FileManager.default.removeItem(at: outURL) } catch { print("Cleanup failed (\(outURL.lastPathComponent)): \(error.localizedDescription)") }
            
            // string
            outURL = sourceDir.appending(path: "Tacitus.txt")
            let _ = try examineString.saveTo(fileName: outURL.path, atomically: true, overwrite: true, separator: ",", encloseElementsBy: nil, asRow: true, stringEncoding: String.Encoding.utf8, compressAsZip: true, zipEntryFileName: nil)
            outURL = sourceDir.appending(path: "Tacitus.txt.zip")
            if let examineString2: SSExamine<String, Double> = try .examine(fromFile: outURL.path, separator: ",",  { token in
                token?.trimmingCharacters(in: .whitespacesAndNewlines)
                
            } )
            {
                print(examineString.isEqual(examineString2))
            }
            else {
                print("error")
            }
            // cleanup Tacitus.txt.zip
            do { try FileManager.default.removeItem(at: outURL) } catch { print("Cleanup failed (\(outURL.lastPathComponent)): \(error.localizedDescription)") }
            
            outURL = sourceDir.appending(path: "Tacitus.json")
            let _ = try examineString.exportJSON(fileName: outURL.path, atomically: true, overwrite: true)
            if let jsonExamine: SSExamine<String, Double> = try .examine(fromJSON: outURL.path) {
                print(examineString.isEqual(jsonExamine))
            }
            // cleanup Tacitus.json
            do { try FileManager.default.removeItem(at: outURL) } catch { print("Cleanup failed (\(outURL.lastPathComponent)): \(error.localizedDescription)") }
        }
    }
    catch let e {
        print(e.localizedDescription)
    }
    
}


private func smokeTests() {
    do {
        
        let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/Zarr.csv", directoryHint: .notDirectory)
        if let examineDoubles: SSExamine<Double, Double> = .loadReals(from: sourceDir.path, separator: ",", encoding: .utf8, allowNaNAndInfinity: true) {
            let zarrdata: [Double] = examineDoubles.itemsAsNumericArray!
            let test1: GoodnessOfFit<Double>.KSTestResult = GoodnessOfFit<Double>.ksTestOneSample(data: zarrdata, testTarget: .normal, boostrapCount: 9999)
            print(test1)
            
            let examineDouble: SSExamine<Double, Double> = .init(usingArray: zarrdata, name: nil, characterSet: nil)
            print(examineDouble.isGaussian!)
            let ca: [Double] = RNGSampler.randoms(n: 150, dist: .weibull(k: 2.3, lambda: 0.5), seed: 0x0312_1965)
            let test: GoodnessOfFit<Double>.KSTestResult = GoodnessOfFit<Double>.ksTestOneSample(data: ca, testTarget: .weibull)
            print(test.pBootstrap)
            let data: [Double] = [10.0, 50.0, 20.0, 5.0, 15.0]
            let lorenz: SSExamine<Double, Double> = .init(usingArray: data, name: nil, characterSet: nil)
            print(lorenz.isGaussian!)
            let res = lorenz.lorenzCurve!
            for point in res {
                print("u = \(point.u), v = \(point.v)")
            }
            
            let charExamine: SSExamine<String, Double> = try! .init(using: chars, levelOfMeasurement: .ordinal, name: nil, characterSet: .alphanumerics)
            print(charExamine.median)
            print(charExamine.computeEntropyProfile())
            let integers: [Int] = [1, 2, 3, 6, 8, 15, 23, 42,33,12,1,233,101,-44]
            let examineInt: SSExamine<Int,Double> = .init(usingArray: integers, name: nil, characterSet: nil)
            print(try examineInt.boxWhisker()!)
            print(examineDouble.Sn!)
            print(examineDouble.Qn!)
            let randoms: Array<Double> = RNGSampler<Double>.randoms(n: 1001, dist: .gaussian(mean: 0, standardDeviation: 1), seed: 0xffe2_c654)
            let randomEx:SSExamine<Double, Double> = SSExamine(usingArray: randoms, name: nil, characterSet: nil)
            print(randomEx.Sn!)
            print(randomEx.Qn!)
            print(randomEx.medianAbsoluteDeviation(center: randomEx.median!)!)
            print(randomEx.meanAbsoluteDeviation(center: randomEx.median!)!)
            print(randomEx.sampleStandardDeviation!)
            print(randomEx.populationStandardDeviation!)
            print("END")
            
            // Example usage of the loader (replaces the previous split.map that caused the error):
            do {
                // Adjust the separator as needed for your file (e.g., ',', ';', ' ', '\t')
                let sourceDir = URL(filePath: #filePath).deletingLastPathComponent().appending(path: "TestFiles/lew.txt", directoryHint: .notDirectory)
                let lewExamine: SSExamine<Double, Double> = .loadReals(from: sourceDir.path, separator: ",")!
                let lewData: [Double] = lewExamine.itemsAsNumericArray!
                print("Loaded \(lewData.count) doubles from file.")
                var times: Array<Double> = []
                var i = 0.0
                for _ in lewData {
                    times.append(i)
                    i = i + 1.0
                }
                let designMatrix = times.map { ti in
                    [1.0, ti]
                }
                let fit = try GLMFit.fit(y: lewData, x: designMatrix)
                
                let dw = TimeSeries<Double,Double>.DurbinWatsonExact.test(residuals: fit.residuals, designMatrix: designMatrix)
                print("Statisitic: \(dw!.statistic)")
                print("p upper tail: \(dw!.upperTail)")
                print("p lower tail: \(dw!.lowerTail)")
                print("p two sided: \(dw!.twoSided)")
                print("ende")
                let gears: [[Double]] = [gear1, gear2, gear3, gear4, gear5, gear6, gear7, gear8, gear9, gear10]
                let allGears: [Double] = gears.flatMap(\.self)
                var examinesGearGrouped:[SSExamine<Double,Double>] = []
                let examineAllGears: SSExamine<Double,Double> = .init(usingArray: allGears, name: "all gears", characterSet: nil)
                i = 1
                let dataFrame: DataFrame<Double, Double> = DataFrame<Double, Double>()
                for a in gears {
                    let e: SSExamine<Double, Double> = .init(usingArray: a, name: "gear" + String(i), characterSet: nil)
                    try dataFrame.append(e, name: e.name)
                    examinesGearGrouped.append(e)
                    i += 1
                }
                if let varTest: Variance<Double>.VarianceEqualityTestResult = try Variance<Double>.bartlettTest(samples: dataFrame, alpha: 0.05) {
                    print(varTest)
                }
                if let varTest: Variance<Double>.VarianceEqualityTestResult = try Variance<Double>.bartlettTest(samples: examinesGearGrouped, alpha: 0.05) {
                    print(varTest)
                }
                if let varTest: Variance<Double>.VarianceEqualityTestResult = try Variance<Double>.leveneTest(samples: dataFrame, alpha: 0.05, testType: .median) {
                    print(varTest)
                }
                if let varTest: Variance<Double>.VarianceEqualityTestResult = try Variance<Double>.leveneTest(samples: examinesGearGrouped, alpha: 0.05, testType: .median) {
                    print(varTest)
                }
                if let varTest: Variance<Double>.ChiSquaredVarianceTestResult = try Variance<Double>.chiSquaredTest(sample: examineAllGears, nominalVariance: 0.01, alpha: 0.05, meanIsKnown: false) {
                    print(varTest)
                }
                if let varTest: Variance<Double>.ChiSquaredVarianceTestResult = try Variance<Double>.chiSquaredTest(sample: allGears, nominalVariance: 0.01, alpha: 0.05, meanIsKnown: false) {
                    print(varTest)
                }
                if let varTest: Variance<Double>.ChiSquaredVarianceTestResult = try Variance<Double>.chiSquaredTest(sample: examineAllGears, nominalVariance: 0.01, alpha: 0.05, meanIsKnown: true) {
                    print(varTest)
                }
                if let varTest: Variance<Double>.ChiSquaredVarianceTestResult = try Variance<Double>.chiSquaredTest(sample: allGears, nominalVariance: 0.01, alpha: 0.05, meanIsKnown: true) {
                    print(varTest)
                }
                for df in dataFrame.data {
                    print( df.isGaussian! )
                }
                let batch1Ex: SSExamine<Double, Double> = .init(usingArray: batch1, name: "ceramic batch1", characterSet: nil)
                let batch2Ex: SSExamine<Double, Double> = .init(usingArray: batch2, name: "ceramic batch2", characterSet: nil)
                if let varTest: Variance<Double>.FRatioTestResult = try Inferential<Double>.HypothesisTesting.VarianceTests.fRatioTest(batch1: batch1Ex, batch2: batch2Ex, alpha: 0.05) {
                    print(varTest)
                }
                if let varTest:  Variance<Double>.FRatioTestResult = try Variance<Double>.fRatioTest(batch1: batch1, batch2: batch2, alpha: 0.05) {
                    print(varTest)
                }
                let quickrun: [Double] = [1,1,2,3,1,1,2,3,5,3,3,2,1,1]
                let quickExamine: SSExamine<Double, Double> = .init(usingArray: quickrun, name: "quick", characterSet: nil)
                if let runsTest: Randomness<Double>.RunsTestResult = try Randomness<Double>.runsTest(data: quickExamine, cuttingPoint: .median, alpha: 0.05, withContinuityCorrection: false) {
                    print(runsTest)
                }
                if let runsTest: Randomness<Double>.RunsTestResult = try Randomness<Double>.runsTest(data: lewExamine, cuttingPoint: .median, alpha: 0.05, withContinuityCorrection: false) {
                    print(runsTest)
                }
                let rr: [Double] = SwiftyStats.RNGSampler<Double>.randoms(n: 1000, dist: .uniform01)
                let rrEx: SSExamine<Double, Double> = .init(usingArray: rr, name: nil, characterSet: nil)
                if let runsTest: Randomness<Double>.RunsTestResult = try Randomness<Double>.runsTest(data: rrEx, cuttingPoint: .median, alpha: 0.05, withContinuityCorrection: false) {
                    print(runsTest)
                }
                var bpv: Double = try Parametric<Double>.binomialTest(successes: 22, numberOfTrials: 100, p0: 0.3, alternative: .less)!
                print(bpv)
                bpv = try Parametric<Double>.binomialTest(successes: 22, numberOfTrials: 100, p0: 0.3, alternative: .greater)!
                print(bpv)
                bpv = try Parametric<Double>.binomialTest(successes: 22, numberOfTrials: 100, p0: 0.3, alternative: .twoSided)!
                print(bpv)
                let binomString = "1i9i2mj3992m3m,29394nm d9394jrn"
                let binomStringExamine: SSExamine<String, Double> = try .init(using: binomString, levelOfMeasurement: .ordinal, name: nil, characterSet: .alphanumerics)
                if let binomTest: Parametric<Double>.BinomialTestResult<String, Double> = try Parametric<Double>.binomialTest(data: charExamine, p0: 0.01, characterSet: .alphanumerics, successCode: "e", alternative: .twoSided) {
                    print(binomTest)
                }
                if let binomTest:  Parametric<Double>.BinomialTestResult<String, Double> = try Parametric<Double>.binomialTest(data: charExamine, p0: 0.01, characterSet: .alphanumerics, successCode: "e", alternative: .less) {
                    print(binomTest)
                }
                if let binomTest:  Parametric<Double>.BinomialTestResult<String, Double> = try Parametric<Double>.binomialTest(data: charExamine, p0: 0.01, characterSet: .alphanumerics, successCode: "e", alternative: .greater) {
                    print(binomTest)
                }
                
                let dataA: [Double] = [0.55, 0.67, 0.43, 0.51, 0.48, 0.60, 0.71, 0.53, 0.44, 0.65, 0.75]
                let dataB: [Double] = [0.49, 0.68, 0.59, 0.72, 0.67, 0.75, 0.65, 0.77, 0.62, 0.48, 0.59]
                print(dataA.sorted())
                print(dataB.sorted())
                for i in 0..<dataA.count {
                    if dataA.contains(dataB[i]) {
                        print(dataB[i])
                    }
                }
                let GA: [String] = ["A","A","A","A","A","A","A","A","A","A","A"]
                let GB: [String] = ["B","B","B","B","B","B","B","B","B","B","B"]
                let cData = dataA + dataB
                let cGroups = GA + GB
                var sorted: GroupedData<Double, String> = .init(data: cData, groups: cGroups)
                print(sorted.sorted())
                var r1: Rank<Double, String, Double> = Rank(sortedData: sorted.sorted().sortedData, sortedGroups: sorted.sorted().sortedGroups)
                print(r1)
                print("A sumOfRanks \(r1.sumOfRanks(g: "A"))")
                print("B sumOfRanks \(r1.sumOfRanks(g: "B"))")
                print("A meanRank \(r1.meanRank(g: "A"))")
                print("B meanRank \(r1.meanRank(g: "B"))")
                print("A sampleSize \(r1.sampleSize(g: "A"))")
                print("B sampleSize \(r1.sampleSize(g: "B"))")
                let exU1: SSExamine<Double, Double> = .init(usingArray: dataA, name: "A", characterSet: nil)
                let exU2: SSExamine<Double, Double> = .init(usingArray: dataB, name: "B", characterSet: nil)
                if let UT = try NonParametric<Double>.mannWhitneyUTest(group1: exU1, group2: exU2, continuityCorrection: true) {
                    print(UT)
                }
                if let UT = try NonParametric<Double>.mannWhitneyUTest(group1: dataA, group2: dataB, continuityCorrection: false) {
                    print(UT)
                }
                let heddA: [Double] = [5,5,8,9,13,13,13,15]
                let heddB: [Double] = [3,3,4,5,5,8,10,16]
                if let UT = try NonParametric<Double>.mannWhitneyUTest(group1: heddA, group2: heddB, continuityCorrection: true) {
                    print(UT)
                }
                if let UT = try NonParametric<Double>.mannWhitneyUTest(group1: heddA, group2: heddB, continuityCorrection: false) {
                    print(UT)
                }
                if let fz = try NonParametric<Double>.powerMannWhitney(effectSize: 2 / 3, alpha: 0.1, alternative: .twoSided, N: 100) {
                    print(fz.n1)
                    print(fz.n2)
                    print(fz.power)
                }
                if let fz = try NonParametric<Double>.nMannWhitney(effectSize: 2 / 3, alpha: 0.1, alternative: .twoSided, power: 0.8) {
                    print(fz.n1)
                    print(fz.n2)
                    print(fz.power)
                }
                // data: Rinne: Taschebuch der Statistik, 4th ed., p581 - sleep data
                let set1: [Double] = [0.7, -1.6, -0.2, -1.2, -0.1, 3.4, 3.7, 0.8, 0.0, 2.0]
                let set2: [Double] = [1.9,  0.8,  1.1,  0.1, -0.1, 4.4, 2.5, 1.6, 4.6, 3.4]
                let exSleep1: SSExamine<Double, Double> = .init(usingArray: set1, name: "sleep1", characterSet: nil)
                let exSleep2: SSExamine<Double, Double> = .init(usingArray: set2, name: "sleep2", characterSet: nil)
                if let st = try NonParametric<Double>.signTest(set1: exSleep1, set2: exSleep2) {
                    print(st)
                }
                let random1: [Double] = RNGSampler<Double>.randoms(n: 1000, dist: .gaussian(mean: 0, standardDeviation: 1))
                let random2: [Double] = RNGSampler<Double>.randoms(n: 1000, dist: .gaussian(mean: 2, standardDeviation: 1))
                if let st = try NonParametric<Double>.signTest(set1: random1, set2: random2) {
                    print(st)
                }
                if let st = try NonParametric<Double>.signTest(data: exSleep1) {
                    print(st)
                }
                if let st = try NonParametric<Double>.signTest(data: exSleep2) {
                    print(st)
                }
                if let p = try Variance<Double>.nFRatioTestBalanced(alpha: 0.05, targetPower: 0.8, theta: 0.9) {
                    print(p)
                }
                // data: Rinne: Taschebuch der Statistik, 4th ed., p552 - steel data
                let steel1: [Double] = [63,65,71,75,90,75,68,74,62,73]
                let steel2: [Double] = [80,78,96,87,88,96,82,73,77,79]
                if let wx = try NonParametric<Double>.wilcoxonSignedRankTest(tier1: steel1, tier2: steel2, zeroHandling: .pratt, continuityCorrection: false, exact: true, exactMaxN: 30) {
                    print(wx)
                }
                let set11: [Double] = [
                    1.83,  0.50,  1.62,  2.48, 1.68, 1.88, 1.55, 3.06, 1.30
                ]
                
                let set21: [Double] = [
                    0.878, 0.647, 0.598, 2.05, 1.06, 1.29, 1.06, 3.14, 1.29
                ]
                if let wx = try NonParametric<Double>.wilcoxonSignedRankTest(tier1: set11, tier2: set21, zeroHandling: .pratt, exact: true, exactMaxN: 30) {
                    print(wx)
                }
                if let ks = try NonParametric<Double>.ksTest(set1: steel1, set2: steel2) {
                    print(ks)
                }
                let x: [Double] = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
                let y: [Double] = [3,4,4,5,6,7,7,8,9,10,10,11,12,13,14,15,15,16,17,18,19]
                if let ww = try NonParametric<Double>.waldWolfowitzTwoSampleTest(set1: x, set2: y, alpha: 0.05, withContinuityCorrection: true) {
                    print(ww)
                }
                let x1: [Double] = [1.2, 0.7, 1.9, 1.1, 0.3, 2.0, 1.4, 0.8, 1.6]
                if let ad = try GoodnessOfFit<Double>.andersonDarlingNormalityTest(data: x1, alpha: 0.05) {
                    print(ad)
                }
                
                let g1:  [Double] = [1.12, 0.95, 1.08, 1.01, 0.99, 1.15, 1.06, 0.92, 1.10, 1.03]
                let g2:  [Double] = [1.20, 1.18, 1.25, 1.14, 1.22, 1.19, 1.28, 1.16, 1.24, 1.21]
                let g3:  [Double] = [1.35, 1.31, 1.29, 1.38, 1.33, 1.36, 1.40, 1.27, 1.34, 1.32]
                let g4:  [Double] = [1.45, 1.50, 1.47, 1.42, 1.55, 1.49, 1.53, 1.44, 1.51, 1.46]
                let g5:  [Double] = [1.60, 1.58, 1.63, 1.66, 1.55, 1.61, 1.64, 1.57, 1.62, 1.59]
                let g6:  [Double] = [1.72, 1.68, 1.75, 1.70, 1.73, 1.67, 1.78, 1.71, 1.74, 1.69]
                let g7:  [Double] = [1.85, 1.80, 1.88, 1.83, 1.90, 1.82, 1.87, 1.81, 1.89, 1.84]
                let g8:  [Double] = [2.00, 1.95, 2.05, 1.98, 2.02, 1.97, 2.08, 1.96, 2.04, 1.99]
                let g9:  [Double] = [2.15, 2.10, 2.18, 2.12, 2.20, 2.11, 2.17, 2.09, 2.19, 2.14]
                let g10: [Double] = [2.30, 2.25, 2.35, 2.28, 2.32, 2.27, 2.38, 2.26, 2.34, 2.29]
                let g11: [Double] = [2.45, 2.40, 2.50, 2.43, 2.48, 2.41, 2.55, 2.42, 2.49, 2.44]
                let gears1: [[Double]] = [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11]
                if let kw = try NonParametric<Double>.kruskalWallisTest(data: gears1, alpha: 0.05) {
                    print(kw)
                }
                let allGears1: [Double] = gears.flatMap(\.self)
                var examinesGearGrouped1:[SSExamine<Double,Double>] = []
                i = 1
                let dataFrame1: DataFrame<Double, Double> = DataFrame<Double, Double>()
                for a in gears1 {
                    let e: SSExamine<Double, Double> = .init(usingArray: a, name: "gear " + String(Int(i)), characterSet: nil)
                    try dataFrame1.append(e, name: e.name)
                    examinesGearGrouped1.append(e)
                    i += 1
                }
                
                if let tk = try PostHoc<Double>.tukeyKramerTest(data: dataFrame1 , alpha: 0.05) {
                    print(PostHoc<Double>.nicePostHocTestResults(results: tk))
                }
                if let tk = try PostHoc<Double>.tukeyKramerTest(data: examinesGearGrouped1, alpha: 0.05) {
                    print(PostHoc<Double>.nicePostHocTestResults(results: tk))
                }
                if let st = try PostHoc<Double>.scheffeTest(data: dataFrame1, alpha: 0.05) {
                    print(PostHoc<Double>.nicePostHocTestResults(results: st))
                }
                
                let gg1 = [6.9,5.4,5.8,4.6,4.0]
                let gg2 = [8.3,6.8,7.8,9.2,6.5]
                let gg3 = [8.0, 10.5,  8.1, 6.9, 9.3]
                let gg4 = [5.8, 3.8 , 6.1 ,5.6 ,6.2 ]
                let allGG = [gg1, gg2, gg3, gg4]
                let dataFrame2: DataFrame<Double, Double> = DataFrame<Double, Double>()
                i = 1
                for d in allGG {
                    let ex = SSExamine<Double,Double>.init(usingArray: d, name: "group 0\(Int(i))", characterSet: nil)
                    try dataFrame2.append(ex, name: ex.name)
                    i += 1
                }
                if let tk = try PostHoc<Double>.tukeyKramerTest(data: dataFrame2 , alpha: 0.05) {
                    print(PostHoc<Double>.nicePostHocTestResults(results: tk))
                }
                if let st = try PostHoc<Double>.scheffeTest(data: dataFrame2, alpha: 0.05) {
                    print(PostHoc<Double>.nicePostHocTestResults(results: st))
                }
                if let an = try Parametric<Double>.oneWayANOVA(data: dataFrame2, postHocTest: .scheffe, alpha: 0.05) {
                    print(an)
                }
                if let an = try Parametric<Double>.oneWayANOVA(data: dataFrame2.data, postHocTest: .tukeyKramer, alpha: 0.05) {
                    print(an)
                }
                let block1: [Double] = [464, 440, 446]
                let block2: [Double] = [441, 393, 334]
                let block3: [Double] = [407, 320, 321]
                let block4: [Double] = [376, 351, 343]
                
                let Fdata: [[Double]] = [block1, block2, block3, block4]
                if let ft = try NonParametric<Double>.friedmanTest(data: Fdata, alpha: 0.05) {
                    print(ft)
                }
                let cars1: [Double] = [18, 15, 18, 16, 17, 15, 14, 14, 14, 15, 15, 14, 15, 14, 22, 18, 21, 21, 10, 10, 11, 9, 28, 25, 19, 16, 17, 19, 18, 14, 14, 14, 14, 12, 13, 13, 18, 22, 19, 18, 23, 26, 25, 20, 21, 13, 14, 15, 14, 17, 11, 13, 12, 13, 15, 13, 13, 14, 22, 28, 13, 14, 13, 14, 15, 12, 13, 13, 14, 13, 12, 13, 18, 16, 18, 18, 23, 11, 12, 13, 12, 18, 21, 19, 21, 15, 16, 15, 11, 20, 21, 19, 15, 26, 25, 16, 16, 18, 16, 13, 14, 14, 14, 28, 19, 18, 15, 15, 16, 15, 16, 14, 17, 16, 15, 18, 21, 20, 13, 23, 20, 23, 18, 19, 25, 26, 18, 16, 16, 15, 22, 22, 24, 23, 29, 25, 20, 18, 19, 18, 27, 13, 17, 13, 13, 13, 30, 26, 18, 17, 16, 15, 18, 21, 19, 19, 16, 16, 16, 16, 25, 26, 31, 34, 36, 20, 19, 20, 19, 21, 20, 25, 21, 19, 21, 21, 19, 18, 19, 18, 18, 18, 30, 31, 23, 24, 22, 20, 22, 20, 21, 17, 18, 17, 18, 17, 16, 19, 19, 36, 27, 23, 24, 34, 35, 28, 29, 27, 34, 32, 28, 26, 24, 19, 28, 24, 27, 27, 26, 24, 30, 39, 35, 34, 30, 22, 27, 20, 18, 28, 27, 34, 31, 29, 27, 24, 23, 38, 36, 25, 38, 26, 22, 36, 27, 27, 32, 28, 31]
                let cars2: [Double] = [24, 27, 27, 25, 31, 35, 24, 19, 28, 23, 27, 20, 22, 18, 20, 31, 32, 31, 32, 24, 26, 29, 24, 24, 33, 33, 32, 28, 19, 32, 34, 26, 30, 22, 22, 33, 39, 36, 28, 27, 21, 24, 30, 34, 32, 38, 37, 30, 31, 37, 32, 47, 41, 45, 34, 33, 24, 32, 39, 35, 32, 37, 38, 34, 34, 32, 33, 32, 25, 24, 37, 31, 36, 36, 34, 38, 32, 38, 32]
                if let tt = try Parametric<Double>.twoSampleTTest(sample1: cars1, sample2: cars2, alpha: 0.05) {
                    print(tt)
                }
                if let tt1 = try Parametric.oneSampleTTest(sample: zarrdata, mean: 5.0, alpha: 0.05) {
                    print(tt1)
                }
                if let mmt = try Parametric.oneWayANOVA(data: examinesGearGrouped1, alpha: 0.05) {
                    print(mmt)
                }
            } catch {
                print(error.localizedDescription)
            }
        }
    }
    catch {
        print(error.localizedDescription)
    }
}
func writeTestValues(examine: SSExamine<Double, Double>, name: String) {
    let sourceDir = URL(filePath: #filePath).deletingLastPathComponent()
    let outURL = sourceDir.appending(path: "TestFiles/RND-\(name).csv", directoryHint: .notDirectory)
    do {
        _ = try examine.saveTo(fileName: outURL.path, atomically: true, overwrite: true, separator: ",", encloseElementsBy: nil, asRow: true, stringEncoding: .utf8, compressAsZip: false, zipEntryFileName: nil)
    }
    catch {
        print(error)
    }

}

