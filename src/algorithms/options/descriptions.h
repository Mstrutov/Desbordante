#pragma once

#include <string>
#include <sstream>

#include "algorithms/metric/enums.h"

namespace algos {

template<typename BetterEnumType>
static std::string EnumToAvailableValues() {
    std::stringstream avail_values;

    avail_values << '[';

    for (auto const& name : BetterEnumType::_names()) {
        avail_values << name << '|';
    }

    avail_values.seekp(-1, avail_values.cur);
    avail_values << ']';

    return avail_values.str();
}

}  // namespace algos

namespace algos::config::descriptions {
constexpr auto kDData = "path to CSV file, relative to ./input_data";
constexpr auto kDSeparator = "CSV separator";
constexpr auto kDHasHeader = "CSV header presence flag [true|false]";
constexpr auto kDEqualNulls = "specify whether two NULLs should be considered equal";
constexpr auto kDThreads = "number of threads to use. If 0, then as many threads are used as the "
                           "hardware can handle concurrently.";
constexpr auto kDError = "error threshold value for Approximate FD algorithms";
constexpr auto kDMaximumLhs = "max considered LHS size";
constexpr auto kDSeed = "RNG seed";
constexpr auto kDMinimumSupport = "minimum support value (between 0 and 1)";
constexpr auto kDMinimumConfidence = "minimum confidence value (between 0 and 1)";
constexpr auto kDInputFormat = "format of the input dataset for AR mining\n[singular|tabular]";
constexpr auto kDTIdColumnIndex = "index of the column where a TID is stored";
constexpr auto kDItemColumnIndex = "index of the column where an item name is stored";
constexpr auto kDFirstColumnTId = "indicates that the first column contains the transaction IDs";
const std::string _kDMetric = "metric to use\n" +
                              EnumToAvailableValues<algos::metric::Metric>();
const auto kDMetric = _kDMetric.c_str();
constexpr auto kDLhsIndices = "LHS column indices for metric FD verification";
constexpr auto kDRhsIndices = "RHS column indices for metric FD verification";
constexpr auto kDParameter = "metric FD parameter";
constexpr auto kDDistFromNullIsInfinity = "specify whether distance from NULL value is infinity "
                                          "(if not, it is 0)";
constexpr auto kDQGramLength = "q-gram length for cosine metric";
const std::string _kDMetricAlgorithm = "MFD algorithm to use\n" +
                                       EnumToAvailableValues<algos::metric::MetricAlgo>();
const auto kDMetricAlgorithm = _kDMetricAlgorithm.c_str();
constexpr auto kDRadius = "maximum difference between a value and the most common value in a "
                          "cluster";
constexpr auto kDRatio = "ratio between the number of deviating values in a cluster and the "
                         "cluster's size";
constexpr auto kDPreciseAlgorithm = "Algorithm that gives exact FDs for typo miner to compare "
                                    "against approximate FDs";
constexpr auto kDApproximateAlgorithm = "Algorithm which gets approximate FDs for typo miner";

constexpr auto kDARAlgorithm = "Association rule mining algorithm implemented by Borgelt";
}
