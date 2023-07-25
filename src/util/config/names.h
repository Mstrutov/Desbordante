#pragma once

namespace util::config::names {
constexpr auto kTable = "table";
constexpr auto kCsvPath = "csv_path";
constexpr auto kSeparator = "separator";
constexpr auto kHasHeader = "has_header";
constexpr auto kEqualNulls = "is_null_equal_null";
constexpr auto kThreads = "threads";
constexpr auto kError = "error";
constexpr auto kMaximumLhs = "max_lhs";
constexpr auto kSeed = "seed";
constexpr auto kMinimumSupport = "minsup";
constexpr auto kMinimumConfidence = "minconf";
constexpr auto kInputFormat = "input_format";
constexpr auto kTIdColumnIndex = "tid_column_index";
constexpr auto kItemColumnIndex = "item_column_index";
constexpr auto kFirstColumnTId = "has_tid";
constexpr auto kMetric = "metric";
constexpr auto kLhsIndices = "lhs_indices";
constexpr auto kRhsIndices = "rhs_indices";
constexpr auto kRhsIndex = "rhs_index";
constexpr auto kParameter = "parameter";
constexpr auto kDistFromNullIsInfinity = "dist_from_null_is_infinity";
constexpr auto kQGramLength = "q";
constexpr auto kMetricAlgorithm = "metric_algorithm";
constexpr auto kRadius = "radius";
constexpr auto kRatio = "ratio";
constexpr auto kPliRelation = "pli_relation";
constexpr auto kBinaryOperation = "bin_operation";
constexpr auto kFuzziness = "fuzziness";
constexpr auto kFuzzinessProbability = "p_fuzz";
constexpr auto kWeight = "weight";
constexpr auto kBumpsLimit = "bumps_limit";
constexpr auto kIterationsLimit = "iterations_limit";
constexpr auto kACSeed = "ac_seed";
constexpr auto kPreciseAlgorithm = "precise_algorithm";
constexpr auto kApproximateAlgorithm = "approximate_algorithm";
constexpr auto kCfdMinimumSupport = "cfd_minsup";
constexpr auto kCfdMinimumConfidence = "cfd_minconf";
constexpr auto kCfdColumnsNumber = "columns_number";
constexpr auto kCfdTuplesNumber = "tuples_number";
constexpr auto kCfdMaximumLhs = "cfd_max_lhs";
constexpr auto kCfdSubstrategy = "cfd_substrategy";
}  // namespace util::config::names
