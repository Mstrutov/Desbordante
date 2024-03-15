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

    static algos::StdParamsMap GetParamMap(CSVConfig const& csv_config) {
        using namespace config::names;
        return {
                {kCsvConfig, csv_config},
                {kError, config::ErrorType{0.0}},
                {kSeed, decltype(algos::pyro::Parameters::seed){0}},
        };
    }

    static algos::StdParamsMap MaxLHSParamMap(tests::TableConfig const& config,
                                              unsigned int max_lhs_) {
        using namespace config::names;
        return {
                {kTable, config.MakeInputTable()},
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
    static std::unique_ptr<algos::FDAlgorithm> CreateAlgorithmInstance(CSVConfig const& config) {
        return algos::CreateAndLoadAlgorithm<T>(GetParamMap(config));
    }

    static std::unique_ptr<algos::FDAlgorithm> CreateMaxLHSAlgorithmInstance(
            tests::TableConfig const& config, unsigned int max_lhs_) {
        return algos::CreateAndLoadAlgorithm<T>(MaxLHSParamMap(config, max_lhs_));
    }

    inline static std::vector<CSVConfigHash> const light_datasets_ = {
            {{tests::kCIPublicHighway10k, 33398},
             {tests::kneighbors10k, 43368},
             {tests::kWDC_astronomical, 22281},
             {tests::kWDC_age, 19620},
             {tests::kWDC_appearances, 25827},
             {tests::kWDC_astrology, 40815},
             {tests::kWDC_game, 6418},
             {tests::kWDC_science, 19620},
             {tests::kWDC_symbols, 28289},
             {tests::kbreast_cancer, 15121},
             {tests::kWDC_kepler, 63730}}};

    inline static std::vector<CSVConfigHash> const heavy_datasets_ = {
            {{tests::kadult, 23075},
             {tests::kCIPublicHighway, 13035},
             {tests::kEpicMeds, 50218},
             {tests::kEpicVitals, 2083},
             {tests::kiowa1kk, 28573},
             {tests::kLegacyPayors, 43612}}};
};
}  // namespace tests
