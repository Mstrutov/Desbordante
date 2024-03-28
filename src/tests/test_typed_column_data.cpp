#include <memory>

#include <gtest/gtest.h>

#include "algorithms/fd/fd_algorithm.h"
#include "all_csv_configs.h"
#include "csv_config_util.h"
#include "model/table/column_layout_typed_relation_data.h"

namespace tests {

namespace mo = model;

using algos::FDAlgorithm;
using mo::TypeId;

struct TypeParsingParams {
    std::vector<mo::TypeId> expected;
    CSVConfig const& csv_config;
};

class TestTypeParsing : public ::testing::TestWithParam<TypeParsingParams> {};

TEST_P(TestTypeParsing, DefaultTest) {
    auto const& [expected, csv_config] = GetParam();
    auto input_table = MakeInputTable(csv_config);
    std::vector<mo::TypedColumnData> column_data{mo::CreateTypedColumnData(*input_table, true)};

    ASSERT_EQ(column_data.size(), expected.size());

    size_t i = 0;
    for (mo::TypedColumnData const& col_data : column_data) {
        EXPECT_EQ(col_data.GetTypeId(), expected[i]) << "Column index: " << i;
        ++i;
    }
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(
    TypeSystem, TestTypeParsing,
    ::testing::Values(
        TypeParsingParams({TypeId::kString, TypeId::kMixed, TypeId::kMixed},
                          kWdcAppearances),
        TypeParsingParams({TypeId::kString, TypeId::kString, TypeId::kString},
                          kWdcAge),
        TypeParsingParams({TypeId::kString, TypeId::kMixed, TypeId::kMixed,
                           TypeId::kMixed},
                          kWdcKepler),
        TypeParsingParams({TypeId::kString, TypeId::kString, TypeId::kMixed,
                           TypeId::kMixed, TypeId::kMixed, TypeId::kString,
                           TypeId::kString, TypeId::kString},
                          kWdcSatellites),
        TypeParsingParams({TypeId::kString, TypeId::kString, TypeId::kInt,
                           TypeId::kInt, TypeId::kInt, TypeId::kInt,
                           TypeId::kInt, TypeId::kUndefined, TypeId::kUndefined,
                           TypeId::kUndefined, TypeId::kUndefined, TypeId::kInt,
                           TypeId::kInt, TypeId::kInt, TypeId::kInt,
                           TypeId::kInt, TypeId::kUndefined, TypeId::kUndefined},
                          kCIPublicHighway700),
        TypeParsingParams({TypeId::kInt, TypeId::kInt, TypeId::kDouble,
                           TypeId::kInt, TypeId::kInt, TypeId::kInt,
                           TypeId::kInt},
                          kNeighbors10k),
        TypeParsingParams({TypeId::kDouble, TypeId::kDouble, TypeId::kDouble,
                           TypeId::kDouble, TypeId::kString},
                          kIris),
        TypeParsingParams({TypeId::kUndefined, TypeId::kUndefined, TypeId::kUndefined,
                           TypeId::kInt, TypeId::kString, TypeId::kDouble,
                           TypeId::kBigInt, TypeId::kMixed, TypeId::kBigInt,
                           TypeId::kMixed, TypeId::kInt},
                          kSimpleTypes),
        TypeParsingParams({TypeId::kString, TypeId::kDate, TypeId::kDate},
                          kACShippingDates)));

// clang-format on

TEST(TypeSystem, SumColumnDoubles) {
    auto input_table = MakeInputTable(kIris);
    std::vector<mo::TypedColumnData> col_data{mo::CreateTypedColumnData(*input_table, true)};
    ASSERT_EQ(col_data.size(), 5);
    mo::TypedColumnData const& col = col_data.front();
    ASSERT_EQ(col.GetTypeId(), static_cast<TypeId>(TypeId::kDouble));
    mo::INumericType const& type = static_cast<mo::INumericType const&>(col.GetType());
    std::unique_ptr<std::byte[]> sum(type.Allocate());
    for (std::byte const* value : col.GetData()) {
        type.Add(sum.get(), value, sum.get());
    }
    mo::Double expected = 876.5;
    EXPECT_DOUBLE_EQ(type.GetValue<mo::Double>(sum.get()), expected);
}

}  // namespace tests
