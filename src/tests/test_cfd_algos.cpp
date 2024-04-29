
// see input_data/cfd_data/LICENSE

#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/cfd/enums.h"
#include "algorithms/cfd/fd_first_algorithm.h"
#include "all_csv_configs.h"
#include "config/names.h"

namespace tests {

static void CheckCfdSetsEquality(std::set<std::string> const& actual,
                                 std::set<std::string> const& expected) {
    ASSERT_EQ(actual.size(), expected.size()) << "count of cfds does not match: expected "
                                              << expected.size() << ", got " << actual.size();

    for (auto const& string_cfd : actual) {
        if (expected.find(string_cfd) == expected.end()) {
            FAIL() << "generated cfd not found in expected";
        }
    }
    SUCCEED();
}

class CFDAlgorithmTest : public ::testing::Test {
protected:
    static std::unique_ptr<algos::cfd::FDFirstAlgorithm> CreateAlgorithmInstance(
            CSVConfig const& csv_config, unsigned minsup, double minconf, char const* substrategy,
            unsigned int max_lhs, size_t columns_number = 0, size_t tuples_number = 0) {
        using namespace config::names;

        algos::StdParamsMap params{
                {kCsvConfig, csv_config},
                {kCfdMinimumSupport, minsup},
                {kCfdMinimumConfidence, minconf},
                {kCfdMaximumLhs, max_lhs},
                {kCfdSubstrategy, algos::cfd::Substrategy::_from_string(substrategy)},
                {kCfdTuplesNumber, tuples_number},
                {kCfdColumnsNumber, columns_number}};
        return algos::CreateAndLoadAlgorithm<algos::cfd::FDFirstAlgorithm>(params);
    }
};

TEST_F(CFDAlgorithmTest, CfdRelationDataStringFormatTest) {
    auto algorithm = CreateAlgorithmInstance(kTennis, 2, 0.85, "dfs", 3, 4, 5);
    algorithm->Execute();
    std::string expected_data =
            "outlook temp humidity windy\nsunny hot high false\nsunny hot high true\n";
    expected_data += "overcast hot high false\nrainy mild high false\nrainy cool normal false\n";
    ASSERT_EQ(algorithm->GetRelationString(), expected_data);
}

TEST_F(CFDAlgorithmTest, CfdRelationDataPartialStringFormatTest) {
    auto algorithm = CreateAlgorithmInstance(kTennis, 8, 0.85, "dfs", 3);
    algorithm->Execute();
    std::vector<int> tids = {0, 2, 4, 6};
    std::string expected_data =
            "outlook temp humidity windy play\nsunny hot high false no\novercast hot high false "
            "yes\n";
    expected_data += "rainy cool normal false yes\novercast cool normal true yes\n";
    ASSERT_EQ(algorithm->GetRelationString(tids), expected_data);
}

TEST_F(CFDAlgorithmTest, FullTennisDataset) {
    auto algorithm = CreateAlgorithmInstance(kTennis, 8, 0.85, "dfs", 3);
    algorithm->Execute();
    std::set<std::string> actual_cfds;
    for (auto const& cfd : algorithm->GetItemsetCfds()) {
        actual_cfds.insert(algorithm->GetCfdString(cfd));
    }
    std::set<std::string> expected_cfds = {"(windy, temp, outlook) => humidity",
                                           "(windy, humidity, outlook) => temp",
                                           "(windy, outlook) => play",
                                           "(outlook, windy=false) => play",
                                           "(windy, temp, outlook) => play",
                                           "(play, temp, outlook) => windy",
                                           "(temp, outlook, play=yes) => windy",
                                           "(play, windy, temp) => outlook",
                                           "(play, temp, windy=false) => outlook",
                                           "(humidity, outlook) => play",
                                           "(humidity, temp, outlook) => play",
                                           "(play, temp, outlook) => humidity",
                                           "(windy, humidity, outlook) => play"};
    CheckCfdSetsEquality(actual_cfds, expected_cfds);

    algorithm = CreateAlgorithmInstance(kTennis, 8, 0.85, "bfs", 3);
    algorithm->Execute();
    CheckCfdSetsEquality(actual_cfds, expected_cfds);
}

TEST_F(CFDAlgorithmTest, PartialMushroomDataset) {
    auto algorithm = CreateAlgorithmInstance(kMushroom, 4, 0.9, "dfs", 4, 4, 50);
    algorithm->Execute();
    std::set<std::string> actual_cfds;
    for (auto const& cfd : algorithm->GetItemsetCfds()) {
        actual_cfds.insert(algorithm->GetCfdString(cfd));
    }
    std::set<std::string> expected_cfds = {"(edible=p) => cap-shape=x",
                                           "(cap-shape=b) => edible=e",
                                           "(cap-color=y) => edible=e",
                                           "(cap-color, edible=p) => cap-shape",
                                           "(edible=p, cap-color=n) => cap-shape=x",
                                           "(cap-surface=f) => edible=e",
                                           "(cap-color, cap-surface=s) => edible",
                                           "(cap-surface, edible=p) => cap-shape",
                                           "(edible=p, cap-surface=y) => cap-shape=x",
                                           "(cap-surface, cap-shape=f) => edible",
                                           "(cap-shape, edible=p, cap-surface=s) => cap-color",
                                           "(cap-color, edible, cap-shape=f) => cap-surface",
                                           "(cap-shape, edible=p, cap-color=w) => cap-surface",
                                           "(edible=p, cap-shape=x, cap-color=w) => cap-surface=y",
                                           "(cap-color, cap-surface, edible=p) => cap-shape",
                                           "(cap-color, cap-surface, cap-shape) => edible",
                                           "(cap-color, cap-shape, cap-surface=s) => edible",
                                           "(cap-color, cap-surface, cap-shape=x) => edible"};

    CheckCfdSetsEquality(actual_cfds, expected_cfds);
}
}  // namespace tests
