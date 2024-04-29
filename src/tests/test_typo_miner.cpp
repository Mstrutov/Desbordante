#include <functional>
#include <utility>

#include <boost/functional/hash.hpp>
#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/pipelines/typo_miner/typo_miner.h"
#include "all_csv_configs.h"
#include "config/names.h"
#include "csv_config_util.h"

namespace tests {
namespace onam = config::names;

struct TestingParam {
    algos::StdParamsMap params;

    TestingParam(CSVConfig const& csv_config, bool const is_null_equal_null, unsigned const max_lhs,
                 double const error, ushort const threads)
        // This dictionary is only created once, so if we put the relation itself here, it
        // won't be reset between different algorithms.
        : params({{onam::kCsvConfig, csv_config},
                  {onam::kEqualNulls, is_null_equal_null},
                  {onam::kMaximumLhs, max_lhs},
                  {onam::kError, error},
                  {onam::kThreads, threads},
                  {onam::kSeed, 0}}) {}
};

/* FD represented as vector where all elements except last are lhs indices and
 * last element is rhs.
 */
using FdByIndices = std::vector<unsigned>;

struct FDsParam : TestingParam {
    std::vector<FdByIndices> expected;

    template <typename... Params>
    explicit FDsParam(std::vector<FdByIndices> expected, Params&&... params)
        : TestingParam(std::forward<Params>(params)...), expected(std::move(expected)) {}
};

struct ClustersParam : TestingParam {
    using ClustersVec = std::vector<model::PLI::Cluster>;
    using FdToClustersMap = std::unordered_map<FdByIndices, ClustersVec, boost::hash<FdByIndices>>;

    FdToClustersMap expected;

    template <typename... Params>
    explicit ClustersParam(FdToClustersMap expected, Params&&... params)
        : TestingParam(std::forward<Params>(params)...), expected(std::move(expected)) {}
};

struct LinesParam : TestingParam {
    using TyposVec = std::vector<model::PLI::Cluster::value_type>;
    using ClusterAndTyposPair = std::pair<model::PLI::Cluster, TyposVec>;
    using FdToTyposMap = std::unordered_map<FdByIndices, std::vector<ClusterAndTyposPair>,
                                            boost::hash<FdByIndices>>;

    double ratio;
    double radius;
    FdToTyposMap expected;

    template <typename... Params>
    explicit LinesParam(FdToTyposMap expected, double ratio, double radius, Params&&... params)
        : TestingParam(std::forward<Params>(params)...),
          ratio(ratio),
          radius(radius),
          expected(std::move(expected)) {}
};

static std::unique_ptr<algos::TypoMiner> ConfToLoadTypoMiner(
        algos::AlgorithmType const algorithm_type, CSVConfig const& csv_config) {
    auto typo_miner = std::make_unique<algos::TypoMiner>(algorithm_type);
    algos::ConfigureFromMap(*typo_miner,
                            algos::StdParamsMap{{onam::kTable, MakeInputTable(csv_config)}});
    return typo_miner;
}

static std::unique_ptr<algos::TypoMiner> CreateTypoMiner(algos::AlgorithmType const algorithm,
                                                         algos::StdParamsMap m) {
    auto typo_miner = std::make_unique<algos::TypoMiner>(algorithm);
    algos::LoadAlgorithm(*typo_miner, m);
    return typo_miner;
}

static std::string MakeJsonFromFds(std::vector<FdByIndices> const& fds) {
    std::string json_fds = "{\"fds\": [";

    for (FdByIndices const& fd : fds) {
        json_fds += "{\"lhs\": [";
        auto last = std::prev(fd.end());
        for (auto it = fd.begin(); it != last; ++it) {
            json_fds += std::to_string(*it);
            if (std::next(it) != last) {
                json_fds += ", ";
            }
        }
        json_fds += "], \"rhs\": ";
        assert(!fd.empty());
        json_fds += std::to_string(fd.back());
        json_fds += "},";
    }

    if (json_fds.back() == ',') {
        json_fds.erase(json_fds.size() - 1);
    }
    json_fds += "]}";

    return json_fds;
}

template <typename F>
static void TestForEachAlgo(F&& test) {
    for (algos::AlgorithmType const algorithm_type :
         algos::GetAllDerived<algos::PliBasedFDAlgorithm>()) {
        if (algorithm_type == +algos::AlgorithmType::pfdtane) continue;
        try {
            test(algorithm_type);
        } catch (std::runtime_error const& e) {
            std::cout << "Exception raised in test: " << e.what() << std::endl;
            FAIL();
        }
    }
}

TEST(SimpleTypoMinerTest, ThrowsOnEmpty) {
    auto const test = [](algos::AlgorithmType const algo) {
        auto typo_miner = ConfToLoadTypoMiner(algo, kTestEmpty);
        ASSERT_THROW(typo_miner->LoadData(), std::runtime_error);
    };
    TestForEachAlgo(test);
}

class ApproxFdsMiningTest : public ::testing::TestWithParam<FDsParam> {};

static void TestEquality(std::string const& expected, std::string const& actual,
                         algos::AlgorithmType const algo) {
    EXPECT_EQ(std::hash<std::string>()(expected), std::hash<std::string>()(actual))
            << "TypoMiner with " << algo._to_string() << " as precise algorithm\n"
            << "Expected:\n\t" << expected << "\nActual:\n\t" << actual << std::endl;
}

TEST_P(ApproxFdsMiningTest, SyntheticTest) {
    auto const test = [&param = GetParam()](algos::AlgorithmType const algo) {
        std::string const expected = MakeJsonFromFds(GetParam().expected);
        auto typo_miner = CreateTypoMiner(algo, GetParam().params);
        typo_miner->Execute();
        std::string actual = typo_miner->GetApproxFDsAsJson();
        TestEquality(expected, actual, algo);
    };
    TestForEachAlgo(test);
}

TEST_P(ApproxFdsMiningTest, ConsistentRepeatedExecution) {
    auto const test = [&param = GetParam()](algos::AlgorithmType const algo) {
        std::string const expected = MakeJsonFromFds(param.expected);
        auto typo_miner = CreateTypoMiner(algo, param.params);
        for (int i = 0; i < 5; ++i) {
            algos::ConfigureFromMap(*typo_miner, algos::StdParamsMap{param.params});
            typo_miner->Execute();
            std::string actual = typo_miner->GetApproxFDsAsJson();
            TestEquality(expected, actual, algo);
        }
    };
    TestForEachAlgo(test);
}

INSTANTIATE_TEST_SUITE_P(TypoMinerTestSuite, ApproxFdsMiningTest,
                         ::testing::Values(
                                 /* Expected fds should be sorted by lhs */
                                 FDsParam({{1, 2}}, kSimpleTypos, true, -1, 0.05, 0),
                                 FDsParam({{0, 1}, {0, 2}, {1, 2}, {0}}, kSimpleTypos, true, -1,
                                          0.81, 0),
                                 FDsParam({{0, 1}, {1, 2}}, kSimpleTypos, true, -1, 0.1, 0)));

class ClustersWithTyposMiningTest : public ::testing::TestWithParam<ClustersParam> {};

static FdByIndices FDtoFdByIndices(FD const& fd) {
    FdByIndices fd_by_indicies = fd.GetLhs().GetColumnIndicesAsVector();
    fd_by_indicies.push_back(fd.GetRhs().GetIndex());
    return fd_by_indicies;
}

TEST_P(ClustersWithTyposMiningTest, FindClustersWithTypos) {
    auto const test = [&param = GetParam()](algos::AlgorithmType const algo) {
        ClustersParam::FdToClustersMap const& expected = GetParam().expected;
        auto typo_miner = CreateTypoMiner(algo, GetParam().params);
        typo_miner->Execute();

        for (FD const& fd : typo_miner->GetApproxFDs()) {
            FdByIndices actual_fd = FDtoFdByIndices(fd);
            ASSERT_TRUE(expected.count(actual_fd));
            std::vector<model::PLI::Cluster> const& expected_cluster = expected.at(actual_fd);
            std::vector<model::PLI::Cluster> const actual_cluster =
                    typo_miner->FindClustersWithTypos(fd);
            EXPECT_EQ(expected_cluster, actual_cluster)
                    << "TyposMiner with " << algo._to_string() << " as precise algorithm\n"
                    << "Clusters of FD: " << fd.ToJSONString();
        }
    };
    TestForEachAlgo(test);
}

INSTANTIATE_TEST_SUITE_P(
        TypoMinerTestSuite, ClustersWithTyposMiningTest,
        ::testing::Values(
                /* Expected clusters should be sorted with respect to TyposMiner::SortCluster
                   sorting. */
                ClustersParam({{FdByIndices{1, 2}, {model::PLI::Cluster{7, 9}}}}, kSimpleTypos,
                              true, -1, 0.05, 0),
                ClustersParam({{FdByIndices{0, 1}, {model::PLI::Cluster{4, 0, 1, 5, 6}}},
                               {FdByIndices{1, 2}, {model::PLI::Cluster{7, 9}}}},
                              kSimpleTypos, true, -1, 0.1, 0),
                ClustersParam({{FdByIndices{0, 1}, {model::PLI::Cluster{4, 0, 1, 5, 6}}},
                               {FdByIndices{1, 2}, {model::PLI::Cluster{7, 9}}},
                               {FdByIndices{0, 2},
                                {model::PLI::Cluster{4, 0, 1, 5, 6}, model::PLI::Cluster{7, 9}}},
                               {FdByIndices{0},
                                {model::PLI::Cluster{8, 2, 3, 7, 9, 0, 1, 4, 5, 6}}}},
                              kSimpleTypos, true, -1, 0.81, 0)));

class SquashClusterTest : public ::testing::TestWithParam<TestingParam> {};

static void VerifySquashed(ColumnLayoutRelationData const& rel, FD const& fd,
                           model::PLI::Cluster const& cluster,
                           std::vector<algos::TypoMiner::SquashedElement> const& squashed) {
    std::vector<int> const& probing_table =
            rel.GetColumnData(fd.GetRhs().GetIndex()).GetProbingTable();
    unsigned cluster_index = 0;
    for (auto const& squashed_element : squashed) {
        int const value = probing_table[squashed_element.tuple_index];
        ASSERT_FALSE(value == model::PLI::kSingletonValueId && squashed_element.amount != 1);
        for (unsigned i = 0; i != squashed_element.amount; ++i, ++cluster_index) {
            /* Check that tuples in one squashed element have equal values in rhs */
            int const actual_value = probing_table[cluster[cluster_index]];
            ASSERT_EQ(value, actual_value)
                    << "Squashed element tuple index: " << squashed_element.tuple_index
                    << "\n\tamount: " << squashed_element.amount << "\nFD: " << fd.ToJSONString()
                    << '\n';
        }
        /* Check that the next tuple (after the last in this squashed element) is not equal to
         * previous one and therefore is correctly not presented at this squashed element */
        ASSERT_TRUE(cluster_index == cluster.size() || value == model::PLI::kSingletonValueId ||
                    value != probing_table[cluster[cluster_index]]);
    }
}

TEST_P(SquashClusterTest, SquashCluster) {
    auto typo_miner = CreateTypoMiner(algos::AlgorithmType::pyro, GetParam().params);
    ColumnLayoutRelationData const& rel = typo_miner->GetRelationData();
    typo_miner->Execute();

    for (FD const& fd : typo_miner->GetApproxFDs()) {
        std::vector<model::PLI::Cluster> const actual_clusters =
                typo_miner->FindClustersWithTypos(fd);
        for (auto const& cluster : actual_clusters) {
            std::vector<algos::TypoMiner::SquashedElement> squashed =
                    typo_miner->SquashCluster(fd, cluster);
            VerifySquashed(rel, fd, cluster, squashed);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(TypoMinerTestSuite, SquashClusterTest,
                         ::testing::Values(TestingParam(kSimpleTypos, true, -1, 0.05, 0),
                                           TestingParam(kSimpleTypos, true, -1, 0.1, 0),
                                           TestingParam(kSimpleTypos, true, -1, 0.81, 0)));

class LinesWithTyposMiningTest : public ::testing::TestWithParam<LinesParam> {};

TEST_P(LinesWithTyposMiningTest, FindLinesWithTypos) {
    auto const& p = GetParam();
    auto typo_miner = CreateTypoMiner(algos::AlgorithmType::pyro, GetParam().params);
    RelationalSchema const* schema = typo_miner->GetRelationData().GetSchema();

    for (auto const& [fd_by_indices, clusters_with_typos] : p.expected) {
        assert(fd_by_indices.size() > 1);
        auto bitset =
                schema->IndicesToBitset(fd_by_indices.cbegin(), std::prev(fd_by_indices.cend()));
        FD fd(schema->GetVertical(std::move(bitset)), schema->GetColumn(fd_by_indices.back()));
        for (auto const& [cluster, typos] : clusters_with_typos) {
            std::vector<model::PLI::Cluster::value_type> const actual =
                    typo_miner->FindLinesWithTypos(fd, cluster, p.radius, p.ratio);
            EXPECT_EQ(typos, actual);
        }
    }
}

INSTANTIATE_TEST_SUITE_P(
        TypoMinerTestSuite, LinesWithTyposMiningTest,
        ::testing::Values(
                LinesParam({{FdByIndices{1, 2},
                             {std::pair(model::PLI::Cluster{7, 9}, std::vector{7, 9})}}},
                           1, -1, kSimpleTypos, true, -1, 0.05, 0),
                LinesParam({{FdByIndices{0, 1},
                             {std::pair(model::PLI::Cluster{4, 0, 1, 5, 6}, std::vector{4})}},
                            {FdByIndices{1, 2},
                             {std::pair(model::PLI::Cluster{7, 9}, std::vector{7, 9})}}},
                           1, -1, kSimpleTypos, true, -1, 0.1, 0)));

}  // namespace tests
