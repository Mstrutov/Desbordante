#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "algo_factory.h"
#include "algorithms/fd/sfd/cords.h"
#include "algorithms/fd/sfd/frequency_handler.h"
#include "algorithms/fd/sfd/sample.h"
#include "all_csv_configs.h"
#include "config/equal_nulls/option.h"
#include "config/max_lhs/type.h"
#include "config/names.h"
#include "csv_config_util.h"
#include "model/table/column_layout_typed_relation_data.h"

namespace {
void AssertVectors(std::vector<Column> const& expected,
                   std::vector<model::ColumnIndex> const& actual) {
    ASSERT_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); i++) {
        ASSERT_EQ(expected[i].GetIndex(), actual[i]);
    }
}

void AssertSFDList(std::list<FD> const& actual,
                   std::vector<std::pair<model::ColumnIndex, model::ColumnIndex>> expected) {
    ASSERT_EQ(actual.size(), expected.size());
    size_t ind = 0;
    for (auto i : actual) {
        ASSERT_EQ(i.GetLhsIndices()[0], expected[ind].first);
        ASSERT_EQ(i.GetRhsIndex(), expected[ind++].second);
    }
}
}  // namespace

namespace tests {

TEST(TestCordsUtils, FrequenciesOfIris) {
    std::vector<std::vector<std::pair<std::string, size_t>>> expected = {
            {{"7.400000", 34}, {"7.600000", 33}, {"4.300000", 32}, {"7.100000", 31},
             {"6.100000", 11}, {"5.400000", 8},  {"6.600000", 25}, {"6.400000", 7},
             {"5.800000", 5},  {"4.800000", 14}, {"5.200000", 18}, {"4.900000", 9},
             {"5.500000", 6},  {"5.700000", 4},  {"4.600000", 16}, {"5.100000", 2},
             {"6.000000", 10}, {"5.600000", 12}, {"4.500000", 27}, {"6.700000", 3},
             {"6.300000", 1},  {"6.500000", 13}, {"6.200000", 15}, {"7.300000", 30},
             {"7.900000", 29}, {"6.900000", 17}, {"6.800000", 20}, {"7.000000", 28},
             {"5.900000", 21}, {"4.700000", 24}, {"5.000000", 0},  {"4.400000", 22},
             {"7.700000", 19}, {"7.200000", 23}, {"5.300000", 26}},

            {{"4.200000", 21}, {"4.400000", 19}, {"4.000000", 18}, {"2.400000", 16},
             {"3.700000", 15}, {"4.100000", 20}, {"3.600000", 14}, {"2.800000", 1},
             {"3.800000", 8},  {"3.200000", 2},  {"3.900000", 17}, {"3.100000", 4},
             {"3.300000", 10}, {"2.900000", 5},  {"2.000000", 22}, {"2.500000", 7},
             {"3.400000", 3},  {"3.500000", 9},  {"2.300000", 12}, {"2.200000", 13},
             {"3.000000", 0},  {"2.700000", 6},  {"2.600000", 11}},

            {{"6.900000", 42}, {"6.300000", 41}, {"6.600000", 40}, {"1.000000", 39},
             {"1.100000", 38}, {"3.000000", 37}, {"3.600000", 34}, {"5.400000", 32},
             {"6.700000", 31}, {"5.300000", 30}, {"5.900000", 29}, {"4.900000", 7},
             {"5.800000", 15}, {"3.700000", 36}, {"5.200000", 27}, {"4.800000", 10},
             {"1.900000", 28}, {"4.500000", 3},  {"1.300000", 4},  {"5.100000", 2},
             {"6.000000", 23}, {"1.600000", 5},  {"4.400000", 12}, {"5.000000", 11},
             {"6.400000", 35}, {"3.500000", 24}, {"1.500000", 0},  {"4.200000", 13},
             {"6.100000", 17}, {"3.800000", 33}, {"1.400000", 1},  {"5.700000", 16},
             {"4.600000", 19}, {"5.500000", 18}, {"1.200000", 22}, {"4.100000", 20},
             {"4.000000", 8},  {"5.600000", 6},  {"3.900000", 21}, {"4.700000", 9},
             {"1.700000", 14}, {"3.300000", 25}, {"4.300000", 26}},

            {{"0.600000", 20}, {"1.700000", 19}, {"1.100000", 18}, {"1.900000", 13},
             {"0.200000", 0},  {"2.400000", 16}, {"1.300000", 1},  {"2.100000", 10},
             {"1.800000", 2},  {"2.200000", 15}, {"0.400000", 6},  {"1.500000", 3},
             {"0.100000", 9},  {"1.400000", 4},  {"2.300000", 5},  {"0.300000", 7},
             {"0.500000", 21}, {"2.500000", 17}, {"1.600000", 14}, {"2.000000", 11},
             {"1.200000", 12}, {"1.000000", 8}},

            {{"Iris-setosa", 2}, {"Iris-versicolor", 1}, {"Iris-virginica", 0}}};

    auto table = MakeInputTable(kIris);
    std::unique_ptr<model::ColumnLayoutTypedRelationData> typed_relation =
            model::ColumnLayoutTypedRelationData::CreateFrom(*table, false);
    std::vector<model::TypedColumnData> const& data = typed_relation->GetColumnData();

    algos::FrequencyHandler handler(data, data.size(), 70);

    ASSERT_EQ(handler.Size(), expected.size());
    for (size_t i = 0; i < expected.size(); i++) {
        ASSERT_EQ(handler.ColumnFrequencyMapSize(i), expected[i].size());
        for (auto const& [value, ordinal_number] : expected[i]) {
            ASSERT_EQ(handler.ContainsValAtColumn(value, i), true);
            ASSERT_EQ(handler.GetValueOrdinalNumberAtColumn(value, i), ordinal_number);
        }
    }
}

TEST(TestCordsUtils, SampleSize) {
    ASSERT_EQ(algos::Sample::CalculateSampleSize(465, 4, 1e-06, 0.05), 4215);
    ASSERT_EQ(algos::Sample::CalculateSampleSize(472, 7, 1e-06, 0.05), 3005);
    ASSERT_EQ(algos::Sample::CalculateSampleSize(4198, 2632, 0.1, 0.1), 111);
    ASSERT_EQ(algos::Sample::CalculateSampleSize(765987, 292784, 0.149804, 0.14367), 168);
    ASSERT_EQ(algos::Sample::CalculateSampleSize(640568, 365981, 0.173673, 0.284799), 64);
    ASSERT_EQ(algos::Sample::CalculateSampleSize(319, 4, 0.0181818, 0.16), 485);
}

class CordsAlgorithmTest : public ::testing::Test {
public:
    struct Config {
        bool is_null_equal_null;
        bool only_sfd;
        long double min_cardinality;
        long double max_diff_vals_proportion;
        long double min_sfd_strength_measure;
        long double min_skew_threshold;
        long double min_structural_zeroes_amount;
        long double max_false_positive_probability;
        long double delta;
        size_t max_amount_of_categories;
        config::MaxLhsType max_lhs;
    };

    static algos::StdParamsMap GetParamMap(CSVConfig const& csv_config, Config const& test_config) {
        using namespace config::names;
        return {{kCsvConfig, csv_config},
                {kEqualNulls, test_config.is_null_equal_null},
                {kOnlySFD, test_config.only_sfd},
                {kMinCard, test_config.min_cardinality},
                {kMaxDiffValsProportion, test_config.max_diff_vals_proportion},
                {kMinSFDStrengthMeasure, test_config.min_sfd_strength_measure},
                {kMinSkewThreshold, test_config.min_skew_threshold},
                {
                        kMinStructuralZeroesAmount,
                        test_config.min_structural_zeroes_amount,
                },
                {kMaxFalsePositiveProbability, test_config.max_false_positive_probability},
                {kDelta, test_config.delta},
                {kMaxAmountOfCategories, test_config.max_amount_of_categories},
                {kMaximumLhs, test_config.max_lhs}};
    }

    static std::unique_ptr<algos::Cords> CreateCordsInstance(CSVConfig const& csv_config,
                                                             Config const& test_config) {
        return algos::CreateAndLoadAlgorithm<algos::Cords>(GetParamMap(csv_config, test_config));
    }
};

using TestConfig = CordsAlgorithmTest::Config;
TestConfig const kTestConfigOnlySfdDefault{
        .is_null_equal_null = true,
        .only_sfd = true,
        .min_cardinality = 0.04L,
        .max_diff_vals_proportion = 0.4L,
        .min_sfd_strength_measure = 0.3L,
        .min_skew_threshold = 0.3L,
        .min_structural_zeroes_amount = 1e-01L,
        .max_false_positive_probability = 1e-06L,
        .delta = 0.05L,
        .max_amount_of_categories = 70LU,
        .max_lhs = 1,
};

TEST_F(CordsAlgorithmTest, LineItem) {
    auto a = CreateCordsInstance(kLineItem, kTestConfigOnlySfdDefault);
    a->Execute();
    AssertVectors(a->GetSoftKeys(), {1, 2, 5, 15});
    AssertVectors(a->GetTrivialColumns(), {});
    AssertSFDList(a->FdList(),
                  {{0, 8},  {0, 9},  {10, 3},  {11, 3},  {12, 3},  {10, 6},  {11, 6},  {12, 6},
                   {10, 7}, {11, 7}, {12, 7},  {8, 9},   {10, 8},  {11, 8},  {12, 8},  {10, 9},
                   {11, 9}, {12, 9}, {10, 13}, {10, 14}, {11, 13}, {11, 14}, {12, 13}, {12, 14}});
}

TEST_F(CordsAlgorithmTest, iris) {
    auto a = CreateCordsInstance(kIris, kTestConfigOnlySfdDefault);

    a->Execute();
    AssertVectors(a->GetSoftKeys(), {});
    AssertVectors(a->GetTrivialColumns(), {});
    AssertSFDList(a->FdList(), {{2, 4}, {3, 4}});
}

TEST_F(CordsAlgorithmTest, CIPublicHighway10k) {
    auto a = CreateCordsInstance(kCIPublicHighway10k, kTestConfigOnlySfdDefault);
    a->Execute();
    AssertVectors(a->GetSoftKeys(), {0, 2});
    AssertVectors(a->GetTrivialColumns(), {1, 7, 8, 9, 10, 14, 15, 16, 17});
    AssertSFDList(a->FdList(), {{11, 3}, {11, 4}, {11, 5}, {11, 6}, {11, 12}});
}

}  // namespace tests
