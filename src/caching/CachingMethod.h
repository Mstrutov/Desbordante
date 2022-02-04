#pragma once

enum class CachingMethod {
    kCoin,
    kNoCaching,
    kAllCaching,
    kEntropy,
    kTrueUniquenessEntropy,
    kMeanEntropyThreshold,
    kHeuristicQ2,
    kGini,
    kInvertedEntropy
};