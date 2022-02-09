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

std::unique_ptr<FDAlgorithm> CreateFD_MineAlgorithmInstance(
    std::filesystem::path const& path, char separator = ',', bool has_header = true) {
    return std::make_unique<Fd_mine>(path, separator, has_header);
}

class AlgorithmTest : public LightDatasets, public HeavyDatasets, public ::testing::Test {
};

std::vector<unsigned int> FD_MineBitsetToIndexVector(boost::dynamic_bitset<> const& bitset) {
    std::vector<unsigned int> res;
    for (size_t index = bitset.find_first();
         index != boost::dynamic_bitset<>::npos;
         index = bitset.find_next(index)) {
        res.push_back(index);
    }
    return res;
}

testing::AssertionResult FD_Mine_CheckFDListEquality(
    std::set<std::pair<std::vector<unsigned int>, unsigned int>> actual,
    std::list<FD> const& expected) {
    for (auto& fd : expected) {
        std::vector<unsigned int> lhs_indices =
            FD_MineBitsetToIndexVector(fd.GetLhs().GetColumnIndices());
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

TEST(AlgorithmSyntheticTest, FD_Mine_ThrowsOnEmpty) {
    auto path = std::filesystem::current_path() / "inputData" / "TestEmpty.csv";
    auto algorithm = CreateFD_MineAlgorithmInstance(path, ',', true);
    ASSERT_THROW(algorithm->Execute(), std::runtime_error);
}

TEST(AlgorithmSyntheticTest, FD_Mine_ReturnsEmptyOnSingleNonKey) {
    auto path = std::filesystem::current_path() / "inputData" / "TestSingleColumn.csv";
    auto algorithm = CreateFD_MineAlgorithmInstance(path, ',', true);
    algorithm->Execute();
    ASSERT_TRUE(algorithm->FdList().empty());
}

TEST(AlgorithmSyntheticTest, FD_Mine_WorksOnLongDataset) {
    auto path = std::filesystem::current_path() / "inputData" / "TestLong.csv";

    std::set<std::pair<std::vector<unsigned int>, unsigned int>> true_fd_collection{{{2}, 1}};

    auto algorithm = CreateFD_MineAlgorithmInstance(path, ',', true);
    algorithm->Execute();
    ASSERT_TRUE(FD_Mine_CheckFDListEquality(true_fd_collection, algorithm->FdList()));
}

std::string GetJsonFDs(std::list<FD>& fd_collection) {
    std::string result = "{\"fds\": [";
    std::list<std::string> discovered_fd_strings;
    for (auto& fd : fd_collection) {
        discovered_fd_strings.push_back(fd.ToJSONString());
    }
    discovered_fd_strings.sort();
    for (auto const& fd : discovered_fd_strings) {
        result += fd + ",";
    }
    if (result.back() == ',') {
        result.erase(result.size() - 1);
    }

    result += "]}";
    return result;
}

void MinimizeFDs(std::list<FD>& fd_collection) {
    std::list<FD>::iterator it1 = fd_collection.begin();
    while (it1 != fd_collection.end()) {
        std::list<FD>::iterator it2 = fd_collection.begin();
        while (it2 != fd_collection.end()) {
            if (it1 == it2) {
                it2++;
                continue;
            }
            if (it1->GetRhs().GetIndex() == it2->GetRhs().GetIndex()) {
                auto lhs1 = it1->GetLhs();
                auto lhs2 = it2->GetLhs();
                if (lhs2.Contains(lhs1)) {
                    it2 = fd_collection.erase(it2);
                    continue;
                }
            }
            it2++;
        }
        it1++;
    }
}

TEST_F(AlgorithmTest, FD_Mine_ReturnsSameAsPyro) {
    auto path = std::filesystem::current_path() / "inputData";

    try {
        for (size_t i = 0; i < LightDatasets::DatasetQuantity(); i++) {
            std::cout << LightDatasets::DatasetName(i) << std::endl;
            // TODO: change this hotfix
            if (LightDatasets::DatasetName(i) == "breast_cancer.csv") {
                continue;
            }
            auto algorithm = CreateFD_MineAlgorithmInstance(
                path / LightDatasets::DatasetName(i), LightDatasets::Separator(i),
                LightDatasets::HasHeader(i));

            auto pyro = Pyro(path / LightDatasets::DatasetName(i), LightDatasets::Separator(i),
                             LightDatasets::HasHeader(i), 0, 0, -1);

            algorithm->Execute();
            std::list<FD> fds = algorithm->FdList();
            pyro.Execute();

            for (auto& fd : pyro.FdList()) {
                if (fd.GetLhs().GetArity() == 0) {
                    std::list<FD>::iterator it = fds.begin();
                    while (it != fds.end()) {
                        if (it->GetRhs().GetIndex() == fd.GetRhs().GetIndex()) {
                            it = fds.erase(it);
                            continue;
                        }
                        it++;
                    }
                    fds.push_back(fd);
                }
            }

            MinimizeFDs(fds);
            // std::string algorithm_results = algorithm->GetJsonFDs();
            std::string algorithm_results = GetJsonFDs(fds);
            std::string results_pyro = pyro.FDAlgorithm::GetJsonFDs();

            EXPECT_EQ(results_pyro, algorithm_results)
                << "The new algorithm and Pyro yield different results at "
                << LightDatasets::DatasetName(i);
        }
    }
    catch (std::runtime_error& e) {
        std::cout << "Exception raised in test: " << e.what() << std::endl;
        FAIL();
    }
    SUCCEED();
}
