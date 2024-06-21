#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/nd/nd.h"
#include "algorithms/nd/nd_verifier/nd_verifier.h"
#include "all_csv_configs.h"
#include "config/names.h"

namespace tests {
namespace onam = config::names;

struct NDVerifyingParams {
    algos::StdParamsMap params;

    NDVerifyingParams(config::IndicesType lhs_indices, config::IndicesType rhs_indices,
                      model::WeightType weight = UINT_MAX, CSVConfig const& csv_config = kTestND,
                      bool null_eq_null = true)
        : params({{onam::kCsvConfig, csv_config},
                  {onam::kLhsIndices, std::move(lhs_indices)},
                  {onam::kRhsIndices, std::move(rhs_indices)},
                  {onam::kEqualNulls, null_eq_null},
                  {onam::kWeight, weight}}) {}
};

class TestNDVerifying : public ::testing::TestWithParam<NDVerifyingParams> {};

TEST_P(TestNDVerifying, DefaultTest) {
    auto const& p = GetParam();
    auto mp = algos::StdParamsMap(p.params);
    auto verifier = algos::CreateAndLoadAlgorithm<algos::nd_verifier::NDVerifier>(mp);
    verifier->Execute();
    EXPECT_EQ(verifier->NDHolds(), true);
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(
        NDVerifierTestSuite, TestNDVerifying,
        ::testing::Values(
            NDVerifyingParams({0}, {1}, 4),
            NDVerifyingParams({0}, {2}, 6),
            NDVerifyingParams({0}, {3}, 4),
            NDVerifyingParams({0}, {4}, 5),
            NDVerifyingParams({0}, {5}, 9),
            NDVerifyingParams({0, 1}, {3, 5}, 3)
            ));

INSTANTIATE_TEST_CASE_P(
        NDVerifierHeavyDatasets, TestNDVerifying,
        ::testing::Values(
            NDVerifyingParams({5}, {6}, 1000000, kIowa1kk),  // I just want to see execution time. Real weight doesn't matter (but it shouldn't be very big)
            NDVerifyingParams({16, 17, 18}, {20, 23}, 1000000, kIowa1kk)  // Also, I want to see how execution time depends on number of columns
        ));

INSTANTIATE_TEST_CASE_P(
    NDVerifierTestNullEqualNull, TestNDVerifying,
    ::testing::Values(  // 6-th column contains 2 values and 7 empty cells
        NDVerifyingParams({0}, {6}, 3, kTestND, true),
        NDVerifyingParams({0}, {6}, 7, kTestND, false)
    ));
// clang-format on

}  // namespace tests
