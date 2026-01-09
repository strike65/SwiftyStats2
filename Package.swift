// swift-tools-version: 6.2
// The swift-tools-version declares the minimum version of Swift required to build this package.

/// Package manifest defining the SwiftyStats library, demo executable, and test targets.

import PackageDescription

let package = Package(
    name: "SwiftyStats",
    platforms: [ .iOS(.v16), .macOS(.v13) ],
    products: [
        .library(
            name: "SwiftyStats",
            targets: ["SwiftyStats"]
        ),
        .executable(
            name: "swiftystats-demo",
            targets: ["SwiftyStatsDemo"]
        ),
    ],
    dependencies: [
        .package(url: "https://github.com/strike65/SwiftyBoost.git", .upToNextMajor(from: "1.0.3")),
        .package(url: "https://github.com/apple/swift-docc-plugin", .upToNextMajor(from: "1.3.0")),
        .package(url: "https://github.com/apple/swift-numerics.git", .upToNextMajor(from: "1.1.1")),
        .package(url: "https://github.com/weichsel/ZIPFoundation.git", .upToNextMajor(from: "0.9.0")),
    ],
    targets: [
        .target(
            name: "SwiftyStatsPrelude",
            dependencies: [
                .product(name: "SwiftyBoost", package: "SwiftyBoost"),
                .product(name: "Numerics", package: "swift-numerics"),
            ],
            path: "Sources/SwiftyStatsPrelude"
        ),
        .target(
            name: "SwiftyStats",
            dependencies: [
                .product(name: "SwiftyBoost", package: "SwiftyBoost"),
                "SwiftyStatsPrelude",
                .product(name: "ZIPFoundation", package: "ZIPFoundation")
            ],
            swiftSettings: [
                .define("ACCELERATE_NEW_LAPACK"),
            ]
        ),
        .executableTarget(
            name: "SwiftyStatsDemo",
            dependencies: ["SwiftyStats"],
            path: "Sources/SwiftyStatsDemo",
            resources: [
                .process("TestFiles")
            ]
        ),
        .testTarget(
            name: "SwiftyStatsTests",
            dependencies: [
                "SwiftyStats",
                .product(name: "SwiftyBoost", package: "SwiftyBoost")
            ],
            path: "Tests/SwiftyStatsTests",
            // By default, sources are expected under Tests/SwiftyStatsTests
        ),
    ]
)
