#include <algorithm>
#include <map>
#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Datasets.h"
#include "Pyro.h"
#include "algorithms/TaneX.h"
#include "RelationalSchema.h"

using ::testing::ContainerEq, ::testing::Eq;

using std::string, std::vector;


/* This is a test suite for algorithm verification. It should be possible to run these tests for any algorithm that:
 * 1. extends FDAlgorithm
 * 2. stores the results in FDAlgorithm::fdCollection_
 * 3. has a constructor with signature: myAlgorithm(fs::path const& path, char separator, bool hasHeader)
 *
 * To test your algorithm, just:
 * 1. include the header
 * 2. in createAlgorithmInstance replace "Tane" with <your algorithm class name>
 * */


std::unique_ptr<FDAlgorithm> createAlgorithmInstance(
        fs::path const& path, char separator = ',', bool hasHeader = true) {
    return std::make_unique<Tane>(path, separator, hasHeader);
}

class AlgorithmTest : public LightDatasets, public HeavyDatasets, public ::testing::Test {
};

std::vector<unsigned int> bitsetToIndexVector(boost::dynamic_bitset<> const& bitset) {
    std::vector<unsigned int> res;
    for (size_t index = bitset.find_first();
         index != boost::dynamic_bitset<>::npos;
         index = bitset.find_next(index)) {
        res.push_back(index);
    }
    return res;
}

testing::AssertionResult checkFDListEquality(
        std::set<std::pair<std::vector<unsigned int>, unsigned int>> actual, std::list<FD> const& expected) {
    for (auto& fd : expected) {
        std::vector<unsigned int> lhsIndices = bitsetToIndexVector(fd.getLhs().getColumnIndices());
        std::sort(lhsIndices.begin(), lhsIndices.end());

        if (auto it = actual.find(std::make_pair(lhsIndices, fd.getRhs().getIndex())); it == actual.end()) {
            return testing::AssertionFailure() << "discovered a false FD: "
                << fd.getLhs().toIndicesString() << "->" << fd.getRhs().toIndicesString();
        } else {
            actual.erase(it);
        }
    }
    return actual.empty() ? testing::AssertionSuccess() : testing::AssertionFailure() << "some FDs remain undiscovered";
}

TEST(AlgorithmSyntheticTest, ReturnsEmptyOnEmpty) {
    auto path = fs::current_path() / "inputData" / "TestEmpty.csv";
    auto algorithm = createAlgorithmInstance(path, ',', true);
    algorithm->execute();
    ASSERT_TRUE(algorithm->fdList().empty());
}

TEST(AlgorithmSyntheticTest, ReturnsEmptyOnSingleNonKey) {
    auto path = fs::current_path() / "inputData" / "TestSingleColumn.csv";
    auto algorithm = createAlgorithmInstance(path, ',', true);
    algorithm->execute();
    ASSERT_TRUE(algorithm->fdList().empty());
}

TEST(AlgorithmSyntheticTest, WorksOnLongDataset) {
    auto path = fs::current_path() / "inputData" / "TestLong.csv";

    std::set<std::pair<std::vector<unsigned int>, unsigned int>> trueFDCollection {{{2}, 1}};

    auto algorithm = createAlgorithmInstance(path, ',', true);
    algorithm->execute();
    ASSERT_TRUE(checkFDListEquality(trueFDCollection, algorithm->fdList()));
}

TEST(AlgorithmSyntheticTest, WorksOnWideDataset) {
    auto path = fs::current_path() / "inputData" / "TestWide.csv";

    std::set<std::pair<std::vector<unsigned int>, unsigned int>> trueFDCollection {
        {{0}, 2},
        {{0}, 4},
        {{2}, 0},
        {{2}, 4},
        {{4}, 0},
        {{4}, 2},
        {{}, 1},
        {{}, 3}
    };

    auto algorithm = createAlgorithmInstance(path, ',', true);
    algorithm->execute();
    ASSERT_TRUE(checkFDListEquality(trueFDCollection, algorithm->fdList()));
}


TEST_F(AlgorithmTest, ReturnsSameFDCollectionHash) {
    auto path = fs::current_path() /"inputData";

    try {
        for (size_t i = 0; i < LightDatasets::datasetQuantity(); i++) {
            auto pyro = Pyro(path / LightDatasets::dataset(i), LightDatasets::separator(i),
                             LightDatasets::hasHeader(i), 0, 0, -1);
            pyro.execute();
            EXPECT_EQ(pyro.FDAlgorithm::fletcher16(), LightDatasets::hash(i))
                                << "FD collection hash changed for " << LightDatasets::dataset(i);
        }
    }
    catch (std::runtime_error& e) {
        std::cout << "Exception raised in test: " << e.what() << std::endl;
        FAIL();
    }
    SUCCEED();
}

TEST_F(AlgorithmTest, ReturnsSameAsPyro) {
    auto path = fs::current_path() /"inputData";

    try {
        for (size_t i = 0; i < LightDatasets::datasetQuantity(); i++) {
            auto algorithm = createAlgorithmInstance(
                    path / LightDatasets::dataset(i),LightDatasets::separator(i),
                    LightDatasets::hasHeader(i));
            algorithm->execute();
            std::string algorithmResults = algorithm->getJsonFDs();
            auto pyro = Pyro(path / LightDatasets::dataset(i), LightDatasets::separator(i),
                             LightDatasets::hasHeader(i), 0, 0, -1);
            pyro.execute();
            std::string resultsPyro = pyro.FDAlgorithm::getJsonFDs();
            // std::cout << "HASH for " << LightDatasets::dataset(i) << ": " << pyro.FDAlgorithm::fletcher16();
            EXPECT_EQ(resultsPyro, algorithmResults)
                << "The new algorithm and Pyro yield different results at " << LightDatasets::dataset(i);
        }
    }
    catch (std::runtime_error& e) {
        std::cout << "Exception raised in test: " << e.what() << std::endl;
        FAIL();
    }
    SUCCEED();
}
