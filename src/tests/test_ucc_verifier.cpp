#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/ucc/ucc_verifier/ucc_verifier.h"
#include "all_csv_configs.h"
#include "config/indices/type.h"
#include "config/names.h"

namespace tests {
namespace onam = config::names;

namespace {

class UCCVerifierSimpleParams {
private:
    algos::StdParamsMap params_map_;
    size_t num_clusters_violating_ucc_ = 0;
    size_t num_rows_violating_ucc_ = 0;
    double expected_error = 0.0;
    std::vector<model::PLI::Cluster> clusters_violating_ucc_;

public:
    UCCVerifierSimpleParams(CSVConfig const& csv_config, config::IndicesType column_indices,
                            size_t const num_clusters_violating_ucc,
                            size_t const num_rows_violating_ucc, double const expected_error,
                            std::vector<model::PLI::Cluster> clusters_violating_ucc)
        : params_map_({{onam::kCsvConfig, csv_config}, {onam::kEqualNulls, true}}),
          num_clusters_violating_ucc_(num_clusters_violating_ucc),
          num_rows_violating_ucc_(num_rows_violating_ucc),
          expected_error(expected_error),
          clusters_violating_ucc_(std::move(clusters_violating_ucc)) {
        // if column_indices is empty then it is not inserted into params_map_, and `ucc_indices`
        // option is not passed to the algorithm. Instead, we assume that the option will be
        // initialized with the default value (set of all table columns)
        if (!column_indices.empty()) {
            params_map_.emplace(onam::kUCCIndices, std::move(column_indices));
        }
    }

    algos::StdParamsMap const& GetParamsMap() const {
        return params_map_;
    }

    size_t GetExpectedNumClustersViolatingUCC() const {
        return num_clusters_violating_ucc_;
    }

    size_t GetExpectedNumRowsViolatingUCC() const {
        return num_rows_violating_ucc_;
    }

    double GetExpectedError() const {
        return expected_error;
    }
    std::vector<model::PLI::Cluster> const& GetExpectedClustersViolatingUCC() const {
        return clusters_violating_ucc_;
    }
};

class TestUCCVerifierSimple : public ::testing::TestWithParam<UCCVerifierSimpleParams> {};

}  // namespace

TEST_P(TestUCCVerifierSimple, DefaultTest) {
    UCCVerifierSimpleParams const& p(GetParam());
    algos::StdParamsMap ucc_verifier_params_map = p.GetParamsMap();
    auto verifier =
            algos::CreateAndLoadAlgorithm<algos::UCCVerifier>(std::move(ucc_verifier_params_map));
    verifier->Execute();
    EXPECT_EQ(verifier->UCCHolds(), p.GetExpectedNumClustersViolatingUCC() == 0);
    EXPECT_EQ(verifier->GetNumRowsViolatingUCC(), p.GetExpectedNumRowsViolatingUCC());
    EXPECT_EQ(verifier->GetNumClustersViolatingUCC(), p.GetExpectedNumClustersViolatingUCC());
    EXPECT_EQ(verifier->GetClustersViolatingUCC(), p.GetExpectedClustersViolatingUCC());
    EXPECT_DOUBLE_EQ(verifier->GetError(), p.GetExpectedError());
}

INSTANTIATE_TEST_SUITE_P(
        UCCVerifierSimpleTestSuite, TestUCCVerifierSimple,
        ::testing::Values(
                UCCVerifierSimpleParams(kTestFD, {0}, 1, 12, 1.0,
                                        {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}}),
                UCCVerifierSimpleParams(kTestFD, {0, 1}, 4, 12, 4.0 * 3.0 * 2.0 / 12.0 / 11.0,
                                        {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}, {9, 10, 11}}),
                UCCVerifierSimpleParams(kTestFD, {0, 1, 2}, 4, 8, 4.0 * 2.0 * 1.0 / 12.0 / 11.0,
                                        {{0, 1}, {3, 4}, {6, 7}, {9, 10}}),
                UCCVerifierSimpleParams(kTestFD, {0, 1, 2, 3, 4, 5}, 3, 6,
                                        3.0 * 2.0 * 1.0 / 12.0 / 11.0, {{3, 4}, {6, 7}, {9, 10}}),
                UCCVerifierSimpleParams(kTestWide, {0}, 0, 0, 0.0, {}),
                UCCVerifierSimpleParams(kTestWide, {0, 1, 2, 3, 4}, 0, 0, 0.0, {}),
                UCCVerifierSimpleParams(kTestFD, {}, 3, 6, 3 * 2.0 * 1.0 / 11.0 / 12.0,
                                        {{3, 4}, {6, 7}, {9, 10}}),
                UCCVerifierSimpleParams(kTestWide, {}, 0, 0, 0.0, {})));

namespace {

class UCCVerifierWithHyUCCParams {
private:
    algos::StdParamsMap ucc_verifier_params_map_;
    algos::StdParamsMap hyucc_params_map_;

public:
    explicit UCCVerifierWithHyUCCParams(CSVConfig const& csv_config)
        : ucc_verifier_params_map_({{onam::kCsvConfig, csv_config}, {onam::kEqualNulls, true}}),
          hyucc_params_map_({{onam::kThreads, static_cast<config::ThreadNumType>(1)},
                             {onam::kCsvConfig, csv_config},
                             {onam::kEqualNulls, true}}) {}

    algos::StdParamsMap const& GetUCCVerifierParamsMap() const {
        return ucc_verifier_params_map_;
    }

    algos::StdParamsMap const& GetHyUCCParamsMap() const {
        return hyucc_params_map_;
    }
};

class TestUCCVerifierWithHyUCC : public ::testing::TestWithParam<UCCVerifierWithHyUCCParams> {};

std::unique_ptr<algos::UCCVerifier> CreateAndExecuteUCCVerifier(
        algos::StdParamsMap params_map, std::vector<unsigned int>&& ucc_candidate) {
    params_map.emplace(onam::kUCCIndices, std::move(ucc_candidate));
    auto verifier = algos::CreateAndLoadAlgorithm<algos::UCCVerifier>(std::move(params_map));
    verifier->Execute();
    return verifier;
}

}  // namespace

TEST_P(TestUCCVerifierWithHyUCC, TestWithHyUCC) {
    UCCVerifierWithHyUCCParams const& p(GetParam());

    algos::StdParamsMap hyucc_params_map = p.GetHyUCCParamsMap();
    auto hyucc = algos::CreateAndLoadAlgorithm<algos::HyUCC>(std::move(hyucc_params_map));
    hyucc->Execute();
    std::list<model::UCC> const& mined_uccs = hyucc->UCCList();

    // run ucc_verifier on each UCC from mined_ucc (UCC must hold)
    for (auto const& current_ucc : mined_uccs) {
        auto verifier = CreateAndExecuteUCCVerifier(p.GetUCCVerifierParamsMap(),
                                                    current_ucc.GetColumnIndicesAsVector());
        EXPECT_TRUE(verifier->UCCHolds());
    }

    // Сase to prevent false positives
    // run ucc_verifier on each UCC from mined_ucc with one column index removed
    // (UCC must not hold because HyUCC returns minimal UCCs)
    for (auto const& current_ucc : mined_uccs) {
        std::vector<unsigned int> current_ucc_vec = current_ucc.GetColumnIndicesAsVector();
        if (current_ucc_vec.size() < 2) {
            continue;
        }
        current_ucc_vec.pop_back();
        auto verifier = CreateAndExecuteUCCVerifier(p.GetUCCVerifierParamsMap(),
                                                    std::move(current_ucc_vec));
        EXPECT_FALSE(verifier->UCCHolds());
    }
}

INSTANTIATE_TEST_SUITE_P(UCCVerifierWithHyUCCTestSuite, TestUCCVerifierWithHyUCC,
                         ::testing::Values(UCCVerifierWithHyUCCParams(kAbalone),
                                           UCCVerifierWithHyUCCParams(kBreastCancer),
                                           UCCVerifierWithHyUCCParams(kCIPublicHighway10k)));
}  // namespace tests
