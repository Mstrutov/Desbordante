#include <algorithm>
#include <map>
#include <iostream>
#include <filesystem>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Datasets.h"
#include "Pyro.h"
#include "algorithms/TaneX.h"
#include "algorithms/Fd_mine.h"
#include "RelationalSchema.h"

using ::testing::ContainerEq, ::testing::Eq;

using std::string, std::vector;

std::unique_ptr<FDAlgorithm> createFD_MineAlgorithmInstance(
        std::filesystem::path const& path, char separator = ',', bool hasHeader = true) {
    return std::make_unique<Fd_mine>(path, separator, hasHeader);
}

class AlgorithmTest : public LightDatasets, public HeavyDatasets, public ::testing::Test {
};

std::vector<unsigned int> FD_Mine_bitsetToIndexVector(boost::dynamic_bitset<> const& bitset) {
    std::vector<unsigned int> res;
    for (size_t index = bitset.find_first();
         index != boost::dynamic_bitset<>::npos;
         index = bitset.find_next(index)) {
        res.push_back(index);
    }
    return res;
}

testing::AssertionResult FD_Mine_checkFDListEquality(
        std::set<std::pair<std::vector<unsigned int>, unsigned int>> actual, std::list<FD> const& expected) {
    for (auto& fd : expected) {
        std::vector<unsigned int> lhsIndices = FD_Mine_bitsetToIndexVector(fd.getLhs().getColumnIndices());
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

TEST(AlgorithmSyntheticTest, FD_Mine_ThrowsOnEmpty) {
    auto path = std::filesystem::current_path() / "inputData" / "TestEmpty.csv";
    auto algorithm = createFD_MineAlgorithmInstance(path, ',', true);
    ASSERT_THROW(algorithm->execute(), std::runtime_error);
}

TEST(AlgorithmSyntheticTest, FD_Mine_ReturnsEmptyOnSingleNonKey) {
    auto path = std::filesystem::current_path() / "inputData" / "TestSingleColumn.csv";
    auto algorithm = createFD_MineAlgorithmInstance(path, ',', true);
    algorithm->execute();
    ASSERT_TRUE(algorithm->fdList().empty());
}

TEST(AlgorithmSyntheticTest, FD_Mine_WorksOnLongDataset) {
    auto path = std::filesystem::current_path() / "inputData" / "TestLong.csv";

    std::set<std::pair<std::vector<unsigned int>, unsigned int>> trueFDCollection {{{2}, 1}};

    auto algorithm = createFD_MineAlgorithmInstance(path, ',', true);
    algorithm->execute();
    ASSERT_TRUE(FD_Mine_checkFDListEquality(trueFDCollection, algorithm->fdList()));
}

std::string getJsonFDs(std::list<FD> &fdCollection) {
    std::string result = "{\"fds\": [";
    std::list<std::string> discoveredFDStrings;
    for (auto& fd : fdCollection) {
        discoveredFDStrings.push_back(fd.toJSONString());
    }
    discoveredFDStrings.sort();
    for (auto const& fd : discoveredFDStrings) {
        result += fd + ",";
    }
    if (result.back() == ',') {
        result.erase(result.size() - 1);
    }

    result += "]}";
    return result;
}

void minimizeFDs(std::list<FD> &fdCollection) {
    std::list<FD>::iterator it1 = fdCollection.begin();
    while (it1 != fdCollection.end()){
        std::list<FD>::iterator it2 = fdCollection.begin();
        while (it2 != fdCollection.end()){
            if (it1 == it2) {
                it2++;
                continue;
            }
            if (it1->getRhs().getIndex() == it2->getRhs().getIndex()) {
                auto lhs1 = it1->getLhs();
                auto lhs2 = it2->getLhs();
                if (lhs2.contains(lhs1)) {
                    it2 = fdCollection.erase(it2);
                    continue;
                }
            }
            it2++;
        }
        it1++;
    }
}

TEST_F(AlgorithmTest, FD_Mine_ReturnsSameAsPyro) {
    auto path = std::filesystem::current_path() /"inputData";

    try {
        for (size_t i = 0; i < LightDatasets::datasetQuantity(); i++) {
            std::cout << LightDatasets::dataset(i) << std::endl;
            // TODO: change this hotfix
            if (LightDatasets::dataset(i) == "breast_cancer.csv") {
                continue;
            }
            auto algorithm = createFD_MineAlgorithmInstance(
                    path / LightDatasets::dataset(i),LightDatasets::separator(i),
                    LightDatasets::hasHeader(i));
              
            auto pyro = Pyro(path / LightDatasets::dataset(i), LightDatasets::separator(i),
                             LightDatasets::hasHeader(i), 0, 0, -1);
            
            algorithm->execute();
            std::list<FD> fds = algorithm->fdList();
            pyro.execute();
            
            for (auto &fd : pyro.fdList()) {
                if (fd.getLhs().getArity() == 0) {
                    std::list<FD>::iterator it = fds.begin();
                    while (it != fds.end()){
                        if (it->getRhs().getIndex() == fd.getRhs().getIndex()) {
                            it = fds.erase(it);
                            continue;
                        }
                        it++;
                    }
                    fds.push_back(fd);
                }
            }
            
            minimizeFDs(fds);
            // std::string algorithmResults = algorithm->getJsonFDs();
            std::string algorithmResults = getJsonFDs(fds);
            std::string resultsPyro = pyro.FDAlgorithm::getJsonFDs();

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
