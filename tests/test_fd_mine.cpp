#include <algorithm>
#include <filesystem>
#include <iostream>
#include <map>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/fd_mine.h"
#include "algorithms/pyro.h"
#include "algorithms/tane.h"
#include "datasets.h"
#include "model/relational_schema.h"
#include "util/config/error/type.h"
#include "util/config/names.h"

using ::testing::ContainerEq, ::testing::Eq;

using algos::FDAlgorithm, algos::Fd_mine, algos::StdParamsMap, util::config::InputTable;

using std::string, std::vector;

namespace onam = util::config::names;

StdParamsMap FD_MineGetParamMap(const std::filesystem::path& path, char separator = ',',
                                bool has_header = true) {
    InputTable parser = std::make_unique<CSVParser>(path, separator, has_header);
    return {{util::config::names::kTable, parser}};
}

std::unique_ptr<FDAlgorithm> ConfToLoadFD_Mine(std::string const& path, char separator = ',',
                                               bool has_header = true) {
    std::unique_ptr<FDAlgorithm> algorithm = std::make_unique<Fd_mine>();
    algos::ConfigureFromMap(*algorithm, FD_MineGetParamMap(path, separator, has_header));
    return algorithm;
}

std::unique_ptr<FDAlgorithm> CreateFD_MineAlgorithmInstance(std::string const& path,
                                                            char separator = ',',
                                                            bool has_header = true) {
    return algos::CreateAndLoadAlgorithm<Fd_mine>(FD_MineGetParamMap(path, separator, has_header));
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

std::set<std::pair<std::vector<unsigned int>, unsigned int>> FD_MineFDsToSet(
        std::list<FD> const& fds) {
    std::set<std::pair<std::vector<unsigned int>, unsigned int>> set;
    for (auto const& fd : fds) {
        auto const& raw_fd = fd.ToRawFD();
        set.emplace(FD_MineBitsetToIndexVector(raw_fd.lhs_), raw_fd.rhs_);
    }
    return set;
}

TEST(AlgorithmSyntheticTest, FD_Mine_ThrowsOnEmpty) {
    auto path = test_data_dir / "TestEmpty.csv";
    auto algorithm = ConfToLoadFD_Mine(test_data_dir / "TestEmpty.csv", ',', true);
    ASSERT_THROW(algorithm->LoadData(), std::runtime_error);
}

TEST(AlgorithmSyntheticTest, FD_Mine_ReturnsEmptyOnSingleNonKey) {
    auto path = test_data_dir / "TestSingleColumn.csv";
    auto algorithm = CreateFD_MineAlgorithmInstance(path, ',', true);
    algorithm->Execute();
    ASSERT_TRUE(algorithm->FdList().empty());
}

TEST(AlgorithmSyntheticTest, FD_Mine_WorksOnLongDataset) {
    auto path = test_data_dir / "TestLong.csv";

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
    namespace onam = util::config::names;

    try {
        for (size_t i = 0; i < LightDatasets::DatasetQuantity(); i++) {
            std::cout << LightDatasets::DatasetName(i) << std::endl;
            // TODO: change this hotfix
            if (LightDatasets::DatasetName(i) == "breast_cancer.csv") {
                continue;
            }
            auto algorithm = CreateFD_MineAlgorithmInstance(
                    test_data_dir / LightDatasets::DatasetName(i), LightDatasets::Separator(i),
                    LightDatasets::HasHeader(i));

            StdParamsMap params_map{{onam::kCsvPath, test_data_dir / LightDatasets::DatasetName(i)},
                                    {onam::kSeparator, LightDatasets::Separator(i)},
                                    {onam::kHasHeader, LightDatasets::HasHeader(i)},
                                    {onam::kSeed, decltype(Parameters::seed){0}},
                                    {onam::kError, util::config::ErrorType{0.0}}};
            auto pyro_ptr = algos::CreateAndLoadAlgorithm<algos::Pyro>(params_map);
            auto& pyro = *pyro_ptr;

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
