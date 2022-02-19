#include <functional>

#include <boost/functional/hash.hpp>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "AlgoFactory.h"
#include "TypoMiner.h"

namespace tests {

struct TestingParam {
    algos::StdParamsMap params;

    TestingParam(std::string const& dataset,
                 char const separator, bool const has_header, bool const is_null_equal_null,
                 unsigned const max_lhs, double const error, ushort const threads)
        : params({{"data", dataset},
                  {"separator", separator},
                  {"has_header", has_header},
                  {"is_null_equal_null", is_null_equal_null},
                  {"max_lhs", max_lhs},
                  {"error", error},
                  {"threads", threads},
                  {"seed", 0}}) {}
};

/* FD represented as vector where all elements except last are lhs indices and
 * last element is rhs.
 */
using FdByIndices = std::vector<unsigned>;

struct FDsParam : TestingParam {
    std::vector<FdByIndices> expected;

    template<typename... Params>
    FDsParam(std::vector<FdByIndices> expected, Params&&... params)
        : TestingParam(std::forward<Params>(params)...), expected(std::move(expected)) {}
};

struct ClustersParam : TestingParam {
    using FdToClustersMap =
        std::unordered_map<FdByIndices, std::vector<util::PLI::Cluster>, boost::hash<FdByIndices>>;
    FdToClustersMap expected;

    template<typename... Params>
    ClustersParam(FdToClustersMap expected, Params&&... params)
        : TestingParam(std::forward<Params>(params)...), expected(std::move(expected)) {}

};

static std::unique_ptr<algos::TypoMiner<FDAlgorithm>> CreateTypoMiner(algos::Algo const algo,
                                                                      algos::StdParamsMap m) {
    auto typo_miner =
        algos::CreateAlgorithmInstance(algos::AlgoMiningType::typos, algo, std::move(m));
    auto casted = static_cast<algos::TypoMiner<FDAlgorithm> *>(typo_miner.release());
    return std::unique_ptr<algos::TypoMiner<FDAlgorithm>>(casted);
}

static std::string MakeJsonFromFds(std::vector<FdByIndices> fds) {
    std::string json_fds = "{\"fds\": [";

    for (FdByIndices const& fd : fds) {
        json_fds += "{lhs: [";
        auto last = std::prev(fd.end());
        for (auto it = fd.begin(); it != last; ++it) {
            json_fds += std::to_string(*it);
            if (std::next(it) != last) {
                json_fds += ", ";
            }
        }
        json_fds += "], rhs: ";
        assert(fd.size() != 0);
        json_fds += std::to_string(fd.back());
        json_fds += "},";
    }

    if (json_fds.back() == ',') {
        json_fds.erase(json_fds.size() - 1);
    }
    json_fds += "]}";

    return json_fds;
}

template<typename F>
static void TestForEachAlgo(F&& test) {
    for (algos::Algo const algo : algos::Algo::_values()) {
        try {
            test(algo);
        } catch (std::runtime_error const& e) {
            std::cout << "Exception raised in test: " << e.what() << std::endl;
            FAIL();
        }
    }
}

TEST(SimpleTypoMinerTest, ThrowsOnEmpty) {
    auto const test = [](algos::Algo const algo) {
        TestingParam p("TestEmpty.csv", ',', true, true, -1, 0.1, 0);
        auto typo_miner = CreateTypoMiner(algo, std::move(p.params));
        ASSERT_THROW(typo_miner->Execute(), std::runtime_error);
    };
    TestForEachAlgo(test);
}

class ApproxFdsMiningTest : public ::testing::TestWithParam<FDsParam> {};

TEST_P(ApproxFdsMiningTest, SyntheticTest) {
    auto const test = [&param = GetParam()](algos::Algo const algo) {
        std::string const expected = MakeJsonFromFds(GetParam().expected);
        auto typo_miner = CreateTypoMiner(algo, GetParam().params);
        typo_miner->Execute();
        std::string actual = typo_miner->GetApproxFDsAsJson();
        EXPECT_EQ(std::hash<std::string>()(expected), std::hash<std::string>()(actual))
            << "TyposMiner with " << algo._to_string() << " as precise algorithm\n"
            << "Expected:\n\t" << expected << "\nActual:\n\t" << actual << std::endl;
    };
    TestForEachAlgo(test);
}

INSTANTIATE_TEST_SUITE_P(
    TypoMinerTestSuite, ApproxFdsMiningTest,
    ::testing::Values(
        /* Expected fds should be sorted by lhs */
        FDsParam({{1, 2}}, "SimpleTypos.csv", ',', true, true, -1, 0.05, 0),
        FDsParam({{0, 1}, {1, 2}}, "SimpleTypos.csv", ',', true, true, -1, 0.1, 0)));


class ClustersWithTyposMiningTest : public ::testing::TestWithParam<ClustersParam> {};

static FdByIndices FDtoFdByIndices(FD const& fd) {
    FdByIndices fd_by_indicies = fd.GetLhs().GetColumnIndicesAsVector();
    fd_by_indicies.push_back(fd.GetRhs().GetIndex());
    return fd_by_indicies;
}

TEST_P(ClustersWithTyposMiningTest, SyntheticTest) {
    auto const test = [&param = GetParam()](algos::Algo const algo) {
        ClustersParam::FdToClustersMap const& expected = GetParam().expected;
        auto typo_miner = CreateTypoMiner(algo, GetParam().params);
        typo_miner->Execute();

        for (FD const& fd : typo_miner->GetApproxFDs()) {
            FdByIndices actual_fd = FDtoFdByIndices(fd);
            ASSERT_TRUE(expected.count(actual_fd));
            std::vector<util::PLI::Cluster> const& expected_cluster = expected.at(actual_fd);
            std::vector<util::PLI::Cluster> const actual_cluster =
                typo_miner->FindClustersWithTypos(fd);
            EXPECT_EQ(expected_cluster, actual_cluster)
                << "TyposMiner with " << algo._to_string() << " as precise algorithm\n"
                << "Clusters of FD: " << fd.ToJSONString() << std::endl;

        }
    };
    TestForEachAlgo(test);
}

INSTANTIATE_TEST_SUITE_P(
    TypoMinerTestSuite, ClustersWithTyposMiningTest,
    ::testing::Values(
        /* Expected clusters should be sorted with respect to TyposMiner::SortCluster sorting. */
        ClustersParam({ {FdByIndices{1, 2}, {util::PLI::Cluster{7, 9}}} },
                      "SimpleTypos.csv", ',', true, true, -1, 0.05, 0),
        ClustersParam({ {FdByIndices{0, 1}, {util::PLI::Cluster{4, 0, 1, 5, 6}}},
                        {FdByIndices{1, 2}, {util::PLI::Cluster{7, 9}}} },
                      "SimpleTypos.csv", ',', true, true, -1, 0.1, 0)));

} // namespace tests
