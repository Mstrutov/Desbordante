#include <gtest/gtest.h>

#include "MetricVerifier.h"
#include "AlgoFactory.h"

namespace tests {

struct MetricVerifyingParams {
    algos::StdParamsMap params;

    MetricVerifyingParams(std::string metric,
                          long double min_parameter,
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
                  {"metric", metric}}) {}
};

class TestMetricVerifying : public ::testing::TestWithParam<MetricVerifyingParams> {};

static std::unique_ptr<algos::MetricVerifier> CreateMetricVerifier(algos::StdParamsMap map) {
    auto metric_verrifier =
        algos::CreateAlgorithmInstance(algos::AlgoMiningType::metric, algos::Algo::metric, std::move(map));
    auto casted = dynamic_cast<algos::MetricVerifier*>(metric_verrifier.release());
    return std::unique_ptr<algos::MetricVerifier>(casted);
}

static bool GetResult(algos::MetricVerifier& metric_verifier) {
    metric_verifier.Execute();
    return metric_verifier.metric_fd_holds;
}

TEST_P(TestMetricVerifying, DefaultTest) {
    std::unique_ptr<algos::MetricVerifier> verifier = CreateMetricVerifier(GetParam().params);

    ASSERT_TRUE(GetResult(*verifier.get()));
    verifier->parameter_ -= verifier->metric_._value == algos::Metric::floating_point
                            ? 1e-9
                            : 1;
    if (verifier->parameter_ < 0) return;
    ASSERT_FALSE(GetResult(*verifier.get()));
}

INSTANTIATE_TEST_SUITE_P(
    MetricVerifierTestSuite, TestMetricVerifying,
    ::testing::Values(
        MetricVerifyingParams("integer", 2, {0, 1}, 2, "TestLong.csv", ',', true),
        MetricVerifyingParams("integer", 1, {0}, 1, "TestLong.csv", ',', true),
        MetricVerifyingParams("integer", 4, {1}, 0, "TestLong.csv", ',', true),
        MetricVerifyingParams("integer", 5, {0}, 2, "TestLong.csv", ',', true),
        MetricVerifyingParams("integer", 0, {2}, 1, "TestLong.csv", ',', true),
        MetricVerifyingParams("integer", 20500, {0}, 4, "TestMetric.csv", ',', true),
        MetricVerifyingParams("integer", 1059, {1}, 4, "TestMetric.csv", ',', true),
        MetricVerifyingParams("integer", 1, {1, 0}, 4, "TestMetric.csv", ',', true),
        MetricVerifyingParams("floating_point", 4.5724231, {0}, 2, "TestMetric.csv", ',', true),
        MetricVerifyingParams("floating_point", 7.53, {0}, 3, "TestMetric.csv", ',', true),
        MetricVerifyingParams("levenshtein", 2, {0}, 5, "TestMetric.csv", ',', true),
        MetricVerifyingParams("levenshtein", 5, {1}, 5, "TestMetric.csv", ',', true),
        MetricVerifyingParams("levenshtein", 4, {0}, 6, "TestMetric.csv", ',', true)
    ));

} // namespace tests