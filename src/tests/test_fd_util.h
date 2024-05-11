#pragma once

#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/fd/fd_algorithm.h"
#include "all_csv_configs.h"
#include "config/error/type.h"
#include "config/names.h"
#include "csv_config_util.h"

namespace tests {
template <typename T>
class AlgorithmTest : public ::testing::Test {
protected:
    static std::unique_ptr<algos::FDAlgorithm> CreateAndConfToLoad(CSVConfig const& csv_config) {
        using namespace config::names;
        using algos::ConfigureFromMap, algos::StdParamsMap;

        std::unique_ptr<algos::FDAlgorithm> algorithm = std::make_unique<T>();
        ConfigureFromMap(*algorithm, StdParamsMap{{kTable, MakeInputTable(csv_config)}});
        return algorithm;
    }

    static algos::StdParamsMap GetParamMap(
            CSVConfig const& csv_config,
            unsigned int max_lhs_ = std::numeric_limits<unsigned int>::max()) {
        using namespace config::names;
        return {
                {kCsvConfig, csv_config},
                {kError, config::ErrorType{0.0}},
                {kSeed, decltype(algos::pyro::Parameters::seed){0}},
                {kMaximumLhs, max_lhs_},
        };
    }

    static void PerformConsistentHashTestOn(std::vector<CSVConfigHash> const& config_hashes) {
        try {
            for (auto const& [csv_config, hash] : config_hashes) {
                auto algorithm = CreateAlgorithmInstance(csv_config);
                algorithm->Execute();
                EXPECT_EQ(algorithm->Fletcher16(), hash)
                        << "FD collection hash changed for " << csv_config.path.filename();
            }
        } catch (std::runtime_error& e) {
            std::cout << "Exception raised in test: " << e.what() << std::endl;
            FAIL();
        }
        SUCCEED();
    }

public:
    static std::unique_ptr<algos::FDAlgorithm> CreateAlgorithmInstance(
            CSVConfig const& config,
            unsigned int max_lhs = std::numeric_limits<unsigned int>::max()) {
        return algos::CreateAndLoadAlgorithm<T>(GetParamMap(config, max_lhs));
    }

    inline static std::vector<CSVConfigHash> const kLightDatasets = {
            {{tests::kCIPublicHighway10k, 33398},
             {tests::kNeighbors10k, 43368},
             {tests::kWdcAstronomical, 22281},
             {tests::kWdcAge, 19620},
             {tests::kWdcAppearances, 25827},
             {tests::kWdcAstrology, 40815},
             {tests::kWdcGame, 6418},
             {tests::kWdcScience, 19620},
             {tests::kWdcSymbols, 28289},
             {tests::kBreastCancer, 15121},
             {tests::kWdcKepler, 63730}}};

    inline static std::vector<CSVConfigHash> const kHeavyDatasets = {
            {{tests::kAdult, 23075},
             {tests::kCIPublicHighway, 13035},
             {tests::kEpicMeds, 50218},
             {tests::kEpicVitals, 2083},
             {tests::kIowa1kk, 28573},
             {tests::kLegacyPayors, 43612}}};
};

template <typename T>
class ApproximateFDTest : public ::testing::Test {
protected:
    // alias for choosing datasets' hashes specialization
    using AlgorithmType = T;

    static std::unique_ptr<algos::FDAlgorithm> CreateAndConfToLoad(CSVConfig const& csv_config) {
        using namespace config::names;
        using algos::ConfigureFromMap, algos::StdParamsMap;

        std::unique_ptr<algos::FDAlgorithm> algorithm = std::make_unique<T>();
        ConfigureFromMap(*algorithm, StdParamsMap{{kTable, MakeInputTable(csv_config)}});
        return algorithm;
    }

    static algos::StdParamsMap GetParamMap(CSVConfig const& csv_config) {
        using namespace config::names;
        return {
             {kCsvConfig, csv_config},
             {kCustomRandom, std::make_pair(true, 47)},
        };
    }

    static void PerformConsistentHashTestOn(std::vector<CSVConfigHash> const& config_hashes) {
        try {
            for (auto const& [csv_config, hash] : config_hashes) {
                auto algorithm = CreateAlgorithmInstance(csv_config);
                algorithm->Execute();
                EXPECT_EQ(algorithm->Fletcher16(), hash)
                        << "FD collection hash changed for " << csv_config.path.filename();
            }
        } catch (std::runtime_error& e) {
            std::cout << "Exception raised in test: " << e.what() << std::endl;
            FAIL();
        }
        SUCCEED();
    }

public:
    static std::unique_ptr<algos::FDAlgorithm> CreateAlgorithmInstance(CSVConfig const& config) {
        return algos::CreateAndLoadAlgorithm<T>(GetParamMap(config));
    }
};

// testing data
template<typename Alg>
struct ApproximateDatasets {
    inline static std::vector<CSVConfigHash> const light_datasets_;
    inline static std::vector<CSVConfigHash> const heavy_datasets_;
};

// specialization fd for EulerFD
template<>
struct ApproximateDatasets<algos::EulerFD> {
    inline static std::vector<CSVConfigHash> const light_datasets_ = {{
             {tests::kCIPublicHighway10k, 33398},
             {tests::kNeighbors10k, 43368},
             {tests::kWdcAstronomical, 2902}, // answer is 9 / 15
             {tests::kWdcAppearances, 64338}, // answer is 3 / 4
             {tests::kWdcAstrology, 40815}, // answer is 34 / 20
             {tests::kWdcSymbols, 28289},
             {tests::kBreastCancer, 15121},
             {tests::kWdcKepler, 17294}, // empty answer, 0 clusters after stripping
    }};

    inline static std::vector<CSVConfigHash> const heavy_datasets_ = {{
             {tests::kAdult, 23075},
             {tests::kCIPublicHighway, 13035},
             {tests::kEpicMeds, 26201}, // answer is 15 / 16
             {tests::kEpicVitals, 2083},
             {tests::kIowa1kk, 57837}, // answer is 2531 / 1584 (average 2k, it is bad seed :(
             {tests::kLegacyPayors, 43612}
    }};
};

}  // namespace tests
