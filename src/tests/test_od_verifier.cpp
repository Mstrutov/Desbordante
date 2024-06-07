#include <gtest/gtest.h>

#include "algo_factory.h"
#include "all_csv_configs.h"
#include "config/names.h"
#include "od/od_verifier/od_verifier.h"

namespace tests {

struct ODVerifyingParams {
    algos::StdParamsMap params;
    size_t const row_violate_ods_by_split = 0;
    size_t const row_violate_ods_by_swap = 0;

    ODVerifyingParams(config::IndicesType lhs_indices, config::IndicesType rhs_indices,
                      config::IndicesType context, bool ascending, size_t const row_error_split = 0,
                      size_t const row_error_swap = 0,
                      CSVConfig const& csv_config = kTestODVerifier)
        : params({{config::names::kCsvConfig, csv_config},
                  {config::names::kLhsIndices, std::move(lhs_indices)},
                  {config::names::kRhsIndices, std::move(rhs_indices)},
                  {config::names::kODContext, std::move(context)},
                  {config::names::kAscendingOD, ascending}}),
          row_violate_ods_by_split(row_error_split),
          row_violate_ods_by_swap(row_error_swap) {}
};

class TestODVerifying : public ::testing::TestWithParam<ODVerifyingParams> {};

TEST_P(TestODVerifying, DefaultTest) {
    auto const& p = GetParam();
    auto mp = algos::StdParamsMap(p.params);
    auto verifier = algos::CreateAndLoadAlgorithm<algos::od_verifier::ODVerifier>(mp);
    verifier->Execute();
    EXPECT_EQ(verifier->GetNumRowsViolateBySwap(), p.row_violate_ods_by_swap);
    EXPECT_EQ(verifier->GetNumRowsViolateBySplit(), p.row_violate_ods_by_split);
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(
        ODVerifierTestSuite, TestODVerifying,
        ::testing::Values(
            ODVerifyingParams({1}, {2}, {0}, true, 0, 0),
            ODVerifyingParams({1}, {2}, {}, true, 1, 2),
            ODVerifyingParams({3}, {4}, {0}, true, 0, 1),
            ODVerifyingParams({1}, {2}, {0}, false, 0, 3),
            ODVerifyingParams({3}, {4}, {0}, false, 0, 2),
            ODVerifyingParams({5}, {6}, {0}, true, 1, 0)
            ));
// clang-format on

}  // namespace tests
