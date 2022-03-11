//#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MetricVerifier.h"
#include "AlgoFactory.h"

struct MetricVerifyingParams {
    algos::StdParamsMap params;

    MetricVerifyingParams(long double min_parameter,
                          std::vector<unsigned int> lhs_indices,
                          unsigned int rhs_index,
                          std::string const& dataset,
                          char const separator,
                          bool const has_header)
        : params({{"parameter", min_parameter},
                  {"lhs_indices", std::move(lhs_indices)},
                  {"rhs_index", rhs_index},
                  {"data", dataset},
                  {"separator", separator},
                  {"has_header", has_header},
                  {"is_null_equal_null", true},
                  {"max_lhs", static_cast<unsigned int>(0)},
                  {"threads", static_cast<ushort>(0)}}) {}
};

class TestMetricVerifying : public ::testing::TestWithParam<MetricVerifyingParams> {};

static std::unique_ptr<MetricVerifier> CreateMetricVerifier(algos::StdParamsMap map) {
    auto metric_verrifier =
        algos::details::CreateMetricVerifier(map);
    auto casted = dynamic_cast<MetricVerifier*>(metric_verrifier.release());
    return std::unique_ptr<MetricVerifier>(casted);
}

static bool GetResult(MetricVerifier& metric_verifier) {
    metric_verifier.Execute();
    return metric_verifier.metric_fd_holds;
}

TEST_P(TestMetricVerifying, DefaultTest) {
    std::unique_ptr<MetricVerifier> verifier = CreateMetricVerifier(GetParam().params);

    ASSERT_TRUE(GetResult(*verifier.get()));
    verifier->parameter_ -= 1e-8;
    if (verifier->parameter_ < 0) return;
    ASSERT_FALSE(GetResult(*verifier.get()));
}

INSTANTIATE_TEST_SUITE_P(
    MetricVerifierTestSuite, TestMetricVerifying,
    ::testing::Values(
        MetricVerifyingParams(2, {0, 1}, 2, "TestLong.csv", ',', true),
        MetricVerifyingParams(1, {0}, 1, "TestLong.csv", ',', true),
        MetricVerifyingParams(4, {1}, 0, "TestLong.csv", ',', true),
        MetricVerifyingParams(5, {0}, 2, "TestLong.csv", ',', true),
        MetricVerifyingParams(0, {2}, 1, "TestLong.csv", ',', true),
        MetricVerifyingParams(20500, {0}, 4, "TestMetric.csv", ',', true),
        MetricVerifyingParams(1059, {1}, 4, "TestMetric.csv", ',', true),
        MetricVerifyingParams(1, {1, 0}, 4, "TestMetric.csv", ',', true),
        MetricVerifyingParams(4.5724231, {0}, 2, "TestMetric.csv", ',', true),
        MetricVerifyingParams(7.53, {0}, 3, "TestMetric.csv", ',', true)
        ));