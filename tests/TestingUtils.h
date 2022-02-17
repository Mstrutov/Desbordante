#pragma once

#include "Datasets.h"
#include "FDAlgorithm.h"

#include <gtest/gtest.h>

#include <filesystem>

template <typename T>
class AlgorithmTest : public LightDatasets, public HeavyDatasets, public ::testing::Test {
protected:
    std::unique_ptr<FDAlgorithm> CreateAlgorithmInstance(
        std::filesystem::path const& path, char separator = ',', bool has_header = true) {
        FDAlgorithm::Config c{ .data = path, .separator = separator, .has_header = has_header };
        c.special_params["error"] = 0.0;
        c.special_params["seed"] = 0;
        return std::make_unique<T>(c);
    }
};
