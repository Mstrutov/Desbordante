#pragma once

#include "Datasets.h"

#include <gtest/gtest.h>

#include <filesystem>

template <typename T>
class AlgorithmTest : public LightDatasets, public HeavyDatasets, public ::testing::Test {
protected:
    std::unique_ptr<FDAlgorithm> CreateAlgorithmInstance(
        std::filesystem::path const& path, char separator = ',', bool has_header = true) {
        return std::make_unique<T>(path, separator, has_header);
    }
};
