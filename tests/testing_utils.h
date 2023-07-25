#pragma once

#include <filesystem>

#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/fd_algorithm.h"
#include "datasets.h"
#include "util/config/error/type.h"
#include "util/config/names.h"

template <typename T>
class AlgorithmTest : public ::testing::Test {
    util::config::InputTable MakeCsvParser(std::string const& path, char separator,
                                           bool has_header) {
        return std::make_shared<CSVParser>(path, separator, has_header);
    }

protected:
    std::unique_ptr<algos::FDAlgorithm> CreateAndConfToLoad(std::string const& path,
                                                            char separator = ',',
                                                            bool has_header = true) {
        using util::config::InputTable, algos::ConfigureFromMap, algos::StdParamsMap;
        std::unique_ptr<algos::FDAlgorithm> algorithm = std::make_unique<T>();
        auto parser = MakeCsvParser(path, separator, has_header);
        ConfigureFromMap(*algorithm, StdParamsMap{{util::config::names::kTable, parser}});
        return algorithm;
    }

    algos::StdParamsMap GetParamMap(const std::filesystem::path& path, char separator = ',',
                                    bool has_header = true) {
        using namespace util::config::names;
        return {
                {kTable, MakeCsvParser(path, separator, has_header)},
                {kError, util::config::ErrorType{0.0}},
                {kSeed, decltype(Parameters::seed){0}},
        };
    }

    std::unique_ptr<algos::FDAlgorithm> CreateAlgorithmInstance(const std::string& filename,
                                                                char separator = ',',
                                                                bool has_header = true) {
        return algos::CreateAndLoadAlgorithm<T>(
                GetParamMap(test_data_dir / filename, separator, has_header));
    }
};
