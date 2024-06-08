#pragma once

#include <sstream>
#include <string>

#include "algorithms/cfd/enums.h"
#include "algorithms/fd/pfdtane/enums.h"
#include "algorithms/metric/enums.h"
#include "util/enum_to_available_values.h"

namespace config::descriptions {

namespace details {
std::string const kDMetricString =
        "metric to use\n" + util::EnumToAvailableValues<algos::metric::Metric>();
std::string const kDMetricAlgorithmString =
        "MFD algorithm to use\n" + util::EnumToAvailableValues<algos::metric::MetricAlgo>();
std::string const kDCfdSubstrategyString = "CFD lattice traversal strategy to use\n" +
                                           util::EnumToAvailableValues<algos::cfd::Substrategy>();
std::string const kDErrorMeasureString =
        "PFD error measure to use\n" + util::EnumToAvailableValues<algos::ErrorMeasure>();
}  // namespace details

constexpr auto kDTable = "table processed by the algorithm";
constexpr auto kDTables = "table collection processed by the algorithm";
constexpr auto kDCsvPath = "path to the CSV table";
constexpr auto kDCsvPaths = "path to the collection of CSV tables";
constexpr auto kDSeparator = "CSV separator";
constexpr auto kDHasHeader = "CSV header presence flag [true|false]";
constexpr auto kDEqualNulls = "specify whether two NULLs should be considered equal";
constexpr auto kDThreads =
        "number of threads to use. If 0, then as many threads are used as the "
        "hardware can handle concurrently.";
constexpr auto kDError = "error threshold value for Approximate FD algorithms";
auto const kDErrorMeasure = details::kDErrorMeasureString.c_str();
constexpr auto kDMaximumLhs = "max considered LHS size";
constexpr auto kDSeed = "RNG seed";
constexpr auto kDMinimumSupport = "minimum support value (between 0 and 1)";
constexpr auto kDMinimumConfidence = "minimum confidence value (between 0 and 1)";
constexpr auto kDInputFormat = "format of the input dataset for AR mining\n[singular|tabular]";
constexpr auto kDTIdColumnIndex = "index of the column where a TID is stored";
constexpr auto kDItemColumnIndex = "index of the column where an item name is stored";
constexpr auto kDFirstColumnTId = "indicates that the first column contains the transaction IDs";
auto const kDMetric = details::kDMetricString.c_str();
constexpr auto kDLhsIndices = "LHS column indices";
constexpr auto kDRhsIndices = "RHS column indices";
constexpr auto kDRhsIndex = "RHS column index";
constexpr auto kDUCCIndices = "column indices for UCC verification";
constexpr auto kDParameter = "metric FD parameter";
constexpr auto kDDistFromNullIsInfinity =
        "specify whether distance from NULL value is infinity "
        "(if not, it is 0)";
constexpr auto kDQGramLength = "q-gram length for cosine metric";
auto const kDMetricAlgorithm = details::kDMetricAlgorithmString.c_str();
constexpr auto kDRadius =
        "maximum difference between a value and the most common value in a "
        "cluster";
constexpr auto kDRatio =
        "ratio between the number of deviating values in a cluster and the "
        "cluster's size";
constexpr auto kDPreciseAlgorithm =
        "Algorithm that gives exact FDs for typo miner to compare "
        "against approximate FDs";
constexpr auto kDApproximateAlgorithm = "Algorithm which gets approximate FDs for typo miner";
constexpr auto kDCfdColumnsNumber =
        "Number of columns in the part of the dataset if you "
        "want to use algo not on the full dataset, but on its part";
constexpr auto kDCfdTuplesNumber =
        "Number of tuples in the part of the dataset if you "
        "want to use algo not on the full dataset, but on its part";
constexpr auto kDCfdMinimumSupport =
        "minimum support value (integer number "
        "between 1 and number of tuples in dataset)";
constexpr auto kDCfdMinimumConfidence = "cfd minimum confidence value (between 0 and 1)";
constexpr auto kDCfdMaximumLhs = "cfd max considered LHS size";
auto const kDCfdSubstrategy = details::kDCfdSubstrategyString.c_str();
constexpr auto kDBinaryOperation = "one of available operations: /, *, +, -";
constexpr auto kDFuzziness = "fraction of exceptional records, lies in (0, 1]";
constexpr auto kDFuzzinessProbability =
        "probability, the fraction of exceptional records that "
        "lie outside the bump intervals is at most Fuzziness, lies in (0, 1]";
constexpr auto kDWeight =
        "value lies in (0, 1]. Closer to 0 - many short intervals. "
        "Closer to 1 - small number of long intervals";
constexpr auto kDBumpsLimit = "max considered intervals amount. Pass 0 to remove limit";
constexpr auto kDTimeLimitSeconds = "max running time of the algorithm. Pass 0 to remove limit";
constexpr auto kDIterationsLimit = "limit for iterations of sampling";
constexpr auto kDACSeed = "seed, needed for choosing a data sample";
constexpr auto kDHllAccuracy = "HyperLogLog approximation accuracy";
constexpr auto kDSampleSize = "Size of a table sample";
constexpr auto kDFindNary = "Detect n-ary inclusion dependencies [true|false]";
constexpr auto kDIgnoreNullCols =
        "Ignore INDs which contain columns filled only with NULLs. May increase "
        "performance but impacts the result. [true|false]";
constexpr auto kDIgnoreConstantCols =
        "Ignore INDs which contain columns filled with only one value. May "
        "increase performance but impacts the result. [true|false]";
constexpr auto kDGraphData = "Path to dot-file with graph";
constexpr auto kDGfdData = "Path to file with GFD";
constexpr auto kDMemLimitMB = "memory limit im MBs";
constexpr auto kDLeftTable = "first table processed by the algorithm";
constexpr auto kDRightTable = "second table processed by the algorithm";
constexpr auto kDPruneNonDisjoint =
        "don't search for dependencies where the LHS decision boundary at the same index as the "
        "RHS decision boundary limits the number of records matched";
constexpr auto kDMinSupport = "minimum support for a dependency's LHS";
constexpr auto kDColumnMatches = "column matches to examine";
constexpr auto kDMaxCardinality = "maximum number of MD matching classifiers";
}  // namespace config::descriptions
