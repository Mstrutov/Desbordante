#include <filesystem>
#include <iostream>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/pyro.h"
#include "datasets.h"
#include "util/config/error/type.h"
#include "util/config/names.h"

namespace tests {

using ::testing::ContainerEq;
using algos::FDAlgorithm, algos::PliBasedFDAlgorithm, algos::StdParamsMap;

namespace fs = std::filesystem;

struct KeysTestParams {
    std::vector<unsigned int> const expected;
    std::string_view const dataset;
    char const sep;
    bool const has_header;
    KeysTestParams(std::vector<unsigned int> expected, std::string_view const dataset,
                   char const sep = ',', bool const has_header = true) noexcept
        : expected(std::move(expected)), dataset(dataset), sep(sep), has_header(has_header) {}
};

class KeysTest : public ::testing::TestWithParam<KeysTestParams> {};

template<typename AlgoInterface>
static inline void GetKeysTestImpl(KeysTestParams const& p) {
    namespace onam = util::config::names;
    auto path = test_data_dir / p.dataset;
    std::vector<unsigned int> actual;
    StdParamsMap params_map{{onam::kCsvPath, path},
                            {onam::kSeparator, p.sep},
                            {onam::kHasHeader, p.has_header},
                            {onam::kSeed, decltype(Parameters::seed){0}},
                            {onam::kError, util::config::ErrorType{0.0}}};
    auto pyro_ptr = algos::CreateAndLoadAlgorithm<algos::Pyro>(params_map);
    auto &pyro = *pyro_ptr;

    pyro.Execute();

    try {
        std::vector<Column const*> const keys = pyro.AlgoInterface::GetKeys();
        actual.reserve(keys.size());
        std::transform(keys.cbegin(), keys.cend(), std::back_inserter(actual),
                       [](Column const* col) { return col->GetIndex(); });
    } catch (std::exception const& e) {
        std::cout << "Exception raised in test: " << e.what() << std::endl;
        FAIL();
    }

    std::sort(actual.begin(), actual.end());

    ASSERT_THAT(actual, ContainerEq(p.expected));
}

TEST_P(KeysTest, FDAlgorithmTest) {
    GetKeysTestImpl<FDAlgorithm>(GetParam());
}

TEST_P(KeysTest, PliBasedAlgorithmTest) {
    GetKeysTestImpl<PliBasedFDAlgorithm>(GetParam());
}

INSTANTIATE_TEST_SUITE_P(, KeysTest,
                         ::testing::Values(KeysTestParams({0, 1, 2}, "WDC_age.csv"),
                                           KeysTestParams({0, 1, 2, 3, 4}, "WDC_game.csv"),
                                           KeysTestParams({0, 2}, "WDC_appearances.csv"),
                                           KeysTestParams({3, 4, 5}, "WDC_astronomical.csv"),
                                           KeysTestParams({0, 2}, "CIPublicHighway700.csv"),
                                           KeysTestParams({}, "abalone.csv", ',', false),
                                           KeysTestParams({}, "adult.csv", ';', false)));
}  // namespace tests
