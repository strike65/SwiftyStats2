//
//  Created by VT on 09.11.25.
//  Copyright © 2025 Volker Thieme. All rights reserved.
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
/// Parameter transforms and constraint handling for maximum likelihood optimization.

import SwiftyStatsPrelude

/// Smooth bijection between unconstrained u-coordinates and constrained theta-parameters.
///
/// Change: Added dToTheta closure which provides the derivative dtheta/du at a given u.
/// This enables a Jacobian-based delta-method covariance transformation from u to theta.
internal struct Transform<T: RealLike>: Sendable {

    let toTheta: @Sendable (T) -> T   // u -> theta
    let toU:     @Sendable (T) -> T   // theta -> u
    let dToTheta: @Sendable (T) -> T  // derivative dtheta/du evaluated at u

    /// Create a transform function pair (and its derivative) for the given parameter constraint.
    static func make(_ c: ParamConstraint<T>) -> Transform<T> {
        switch c {
        case .real:
            // theta = u, dtheta/du = 1
            return .init(
                toTheta: { $0 },
                toU: { $0 },
                dToTheta: { _ in T.one }
            )
        case .positive:
            // theta = softplus(u), u = softplus^{-1}(theta)
            // dtheta/du = logistic(u)
            return .init(
                toTheta: { Helpers.softplus($0) },
                toU: { th in
                    precondition(th > 0, "theta must be > 0")
                    return Helpers.invSoftplus(th)
                },
                dToTheta: { u in Helpers.logistic(u) }
            )
        case .unitInterval:
            // theta = logistic(u), u = logit(theta)
            // dtheta/du = logistic(u) * (1 - logistic(u))
            return .init(
                toTheta: { Helpers.logistic($0) },
                toU: { th in
                    precondition(th > 0 && th < 1, "theta in (0, 1) required")
                    return Helpers.logit(th)
                },
                dToTheta: { u in
                    let s = Helpers.logistic(u)
                    return s * (T.one - s)
                }
            )
        case .lowerBound(let a):
            // theta = a + softplus(u), u = softplus^{-1}(theta - a), dtheta/du = logistic(u)
            return .init(
                toTheta: { a + Helpers.softplus($0) },
                toU: { th in
                    precondition(th > a, "theta > a required")
                    return Helpers.invSoftplus(th - a)
                },
                dToTheta: { u in Helpers.logistic(u) }
            )
        case .upperBound(let b):
            // theta = b - softplus(u), u = softplus^{-1}(b - theta), dtheta/du = - logistic(u)
            return .init(
                toTheta: { b - Helpers.softplus($0) },
                toU: { th in
                    precondition(th < b, "theta < b required")
                    return Helpers.invSoftplus(b - th)
                },
                dToTheta: { u in -Helpers.logistic(u) }
            )
        case .interval(let a, let b):
            precondition(a < b, "interval requires a < b")
            // theta = a + (b-a) * logistic(u)
            // dtheta/du = (b-a) * logistic(u) * (1 - logistic(u))
            return .init(
                toTheta: { a + (b - a) * Helpers.logistic($0) },
                toU: { th in
                    precondition(th > a && th < b, "theta in (a, b) required")
                    let t = (th - a) / (b - a)
                    return Helpers.logit(t)
                },
                dToTheta: { u in
                    let s = Helpers.logistic(u)
                    return (b - a) * s * (T.one - s)
                }
            )
        }
    }
}

