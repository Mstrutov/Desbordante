#include <algorithm>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "algorithms/fd/depminer/depminer.h"
#include "algorithms/fd/dfd/dfd.h"
#include "algorithms/fd/fastfds/fastfds.h"
#include "algorithms/fd/fdep/fdep.h"
#include "algorithms/fd/fun/fun.h"
#include "algorithms/fd/hyfd/hyfd.h"
#include "algorithms/fd/pyro/pyro.h"
#include "algorithms/fd/tane/tane.h"
#include "model/table/relational_schema.h"
#include "test_fd_util.h"

using std::string, std::vector;
using ::testing::ContainerEq, ::testing::Eq;

namespace fs = std::filesystem;

namespace tests {

/* This is a test suite for algorithm verification. It should be possible to run these tests for any
 * algorithm that:
 * 1. extends FDAlgorithm
 * 2. stores the results in FDAlgorithm::fd_collection_
 * 3. has a constructor with signature: myAlgorithm(fs::path const& path, char separator, bool
 * has_header)
 *
 * To test your algorithm, just:
 * 1. include the header
 * 2. in CreateAlgorithmInstance replace "Tane" with <your algorithm class name>
 * */

std::vector<unsigned int> BitsetToIndexVector(boost::dynamic_bitset<> const& bitset) {
    std::vector<unsigned int> res;
    for (size_t index = bitset.find_first(); index != boost::dynamic_bitset<>::npos;
         index = bitset.find_next(index)) {
        res.push_back(index);
    }
    return res;
}

testing::AssertionResult CheckFdListEquality(
        std::set<std::pair<std::vector<unsigned int>, unsigned int>> actual,
        std::list<FD> const& expected) {
    for (auto& fd : expected) {
        std::vector<unsigned int> lhs_indices = BitsetToIndexVector(fd.GetLhs().GetColumnIndices());
        std::sort(lhs_indices.begin(), lhs_indices.end());

        if (auto it = actual.find(std::make_pair(lhs_indices, fd.GetRhs().GetIndex()));
            it == actual.end()) {
            return testing::AssertionFailure()
                   << "discovered a false FD: " << fd.GetLhs().ToIndicesString() << "->"
                   << fd.GetRhs().ToIndicesString();
        } else {
            actual.erase(it);
        }
    }
    return actual.empty() ? testing::AssertionSuccess()
                          : testing::AssertionFailure() << "some FDs remain undiscovered";
}

std::set<std::pair<std::vector<unsigned int>, unsigned int>> FDsToSet(std::list<FD> const& fds) {
    std::set<std::pair<std::vector<unsigned int>, unsigned int>> set;
    for (auto const& fd : fds) {
        auto const& raw_fd = fd.ToRawFD();
        set.emplace(BitsetToIndexVector(raw_fd.lhs_), raw_fd.rhs_);
    }
    return set;
}

TYPED_TEST_SUITE_P(AlgorithmTest);

TYPED_TEST_P(AlgorithmTest, ThrowsOnEmpty) {
    auto algorithm = TestFixture::CreateAndConfToLoad(kTestEmpty);
    ASSERT_THROW(algorithm->LoadData(), std::runtime_error);
}

TYPED_TEST_P(AlgorithmTest, ReturnsEmptyOnSingleNonKey) {
    auto algorithm = TestFixture::CreateAlgorithmInstance(kTestSingleColumn);
    algorithm->Execute();
    ASSERT_TRUE(algorithm->FdList().empty());
}

TYPED_TEST_P(AlgorithmTest, WorksOnLongDataset) {
    std::set<std::pair<std::vector<unsigned int>, unsigned int>> true_fd_collection{{{2}, 1}};

    auto algorithm = TestFixture::CreateAlgorithmInstance(kTestLong);
    algorithm->Execute();
    ASSERT_TRUE(CheckFdListEquality(true_fd_collection, algorithm->FdList()));
}

TYPED_TEST_P(AlgorithmTest, WorksOnWideDataset) {
    std::set<std::pair<std::vector<unsigned int>, unsigned int>> true_fd_collection{
            {{0}, 2}, {{0}, 4}, {{2}, 0}, {{2}, 4}, {{4}, 0}, {{4}, 2}, {{}, 1}, {{}, 3}};

    auto algorithm = TestFixture::CreateAlgorithmInstance(kTestWide);
    algorithm->Execute();
    ASSERT_TRUE(CheckFdListEquality(true_fd_collection, algorithm->FdList()));
}

TYPED_TEST_P(AlgorithmTest, LightDatasetsConsistentHash) {
    TestFixture::PerformConsistentHashTestOn(TestFixture::kLightDatasets);
}

TYPED_TEST_P(AlgorithmTest, HeavyDatasetsConsistentHash) {
    TestFixture::PerformConsistentHashTestOn(TestFixture::kHeavyDatasets);
}

TYPED_TEST_P(AlgorithmTest, ConsistentRepeatedExecution) {
    auto algorithm = TestFixture::CreateAlgorithmInstance(kWdcAstronomical);
    algorithm->Execute();
    auto first_res = FDsToSet(algorithm->FdList());
    for (int i = 0; i < 3; ++i) {
        algos::ConfigureFromMap(*algorithm, TestFixture::GetParamMap(kWdcAstronomical));
        algorithm->Execute();
        ASSERT_TRUE(CheckFdListEquality(first_res, algorithm->FdList()));
    }
}

namespace {
void MaxLhsTestFun(CSVConfig config, std::list<FD> const& fds_list, config::MaxLhsType max_lhs) {
    using namespace config::names;
    algos::StdParamsMap verify_params = {
            {kCsvConfig, config},
            {kError, config::ErrorType{0.0}},
            {kMaximumLhs, max_lhs},
    };
    auto verify_algo = algos::CreateAndLoadAlgorithm<algos::Pyro>(verify_params);
    verify_algo->Execute();
    auto verify_list = FDsToSet(verify_algo->FdList());
    ASSERT_TRUE(CheckFdListEquality(verify_list, fds_list));
    for (auto& fd : fds_list) {
        ASSERT_TRUE(fd.GetLhs().GetArity() <= max_lhs);
    }
}
}  // namespace

TYPED_TEST_P(AlgorithmTest, MaxLHSOptionWork) {
    config::MaxLhsType max_lhs = 2;

    auto algo = TestFixture::CreateAlgorithmInstance(kTestFD, max_lhs);
    algo->Execute();
    MaxLhsTestFun(kTestFD, algo->FdList(), max_lhs);

    auto algo_large = TestFixture::CreateAlgorithmInstance(kCIPublicHighway700, max_lhs);
    algo_large->Execute();
    MaxLhsTestFun(kCIPublicHighway700, algo_large->FdList(), max_lhs);
}

REGISTER_TYPED_TEST_SUITE_P(AlgorithmTest, ThrowsOnEmpty, ReturnsEmptyOnSingleNonKey,
                            WorksOnLongDataset, WorksOnWideDataset, LightDatasetsConsistentHash,
                            HeavyDatasetsConsistentHash, ConsistentRepeatedExecution,
                            MaxLHSOptionWork);

using Algorithms =
        ::testing::Types<algos::Tane, algos::Pyro, algos::FastFDs, algos::DFD, algos::Depminer,
                         algos::FDep, algos::FUN, algos::hyfd::HyFD, algos::PFDTane>;
INSTANTIATE_TYPED_TEST_SUITE_P(AlgorithmTest, AlgorithmTest, Algorithms);

}  // namespace tests
