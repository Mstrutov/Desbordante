#include <filesystem>
#include <memory>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ColumnLayoutTypedRelationData.h"
#include "CSVParser.h"

namespace tests {

namespace fs = std::filesystem;
namespace mo = model;

using mo::TypeId;

struct TypeParsingParams {
    std::vector<mo::TypeId> expected;
    std::string_view dataset;
    char const sep;
    bool const has_header;
    TypeParsingParams(std::vector<mo::TypeId> expected,
                      std::string_view dataset,
                      char const sep = ',',
                      bool const has_header = true) noexcept
        : expected(std::move(expected)), dataset(dataset),
          sep(sep), has_header(has_header) {}
};

class TestTypeParsing : public ::testing::TestWithParam<TypeParsingParams> {};

static inline std::vector<mo::TypedColumnData> CreateColumnData(std::string_view data, char sep,
                                                                bool has_header) {
    auto const path = fs::current_path() / "inputData" / data;
    CSVParser input_generator(path, sep, has_header);
    std::unique_ptr<model::ColumnLayoutTypedRelationData> relation_data =
        model::ColumnLayoutTypedRelationData::CreateFrom(input_generator, true, -1, -1);
    std::vector<mo::TypedColumnData> col_data = std::move(relation_data->GetColumnData());
    return col_data;
}

TEST_P(TestTypeParsing, DefaultTest) {
    TypeParsingParams const& p = GetParam();
    std::vector<mo::TypedColumnData> column_data = CreateColumnData(p.dataset, p.sep, p.has_header);

    ASSERT_EQ(column_data.size(), p.expected.size());

    size_t i = 0;
    for (mo::TypedColumnData const& col_data : column_data) {
        EXPECT_EQ(col_data.GetTypeId(), p.expected[i]) << "Column index: " << i;
        ++i;
    }
}

INSTANTIATE_TEST_SUITE_P(
    TypeSystem, TestTypeParsing,
    ::testing::Values(
        TypeParsingParams({TypeId::kString, TypeId::kMixed, TypeId::kMixed},
                          "WDC_appearances.csv"),
        TypeParsingParams({TypeId::kString, TypeId::kString, TypeId::kString},
                          "WDC_age.csv"),
        TypeParsingParams({TypeId::kString, TypeId::kMixed, TypeId::kMixed,
                           TypeId::kMixed},
                          "WDC_kepler.csv"),
        TypeParsingParams({TypeId::kString, TypeId::kString, TypeId::kMixed,
                           TypeId::kMixed, TypeId::kMixed, TypeId::kString,
                           TypeId::kString, TypeId::kString},
                          "WDC_satellites.csv"),
        TypeParsingParams({TypeId::kString, TypeId::kString, TypeId::kInt,
                           TypeId::kInt, TypeId::kInt, TypeId::kInt,
                           TypeId::kInt, TypeId::kUndefined, TypeId::kUndefined,
                           TypeId::kUndefined, TypeId::kUndefined, TypeId::kInt,
                           TypeId::kInt, TypeId::kInt, TypeId::kInt,
                           TypeId::kInt, TypeId::kUndefined, TypeId::kUndefined},
                          "CIPublicHighway700.csv"),
        TypeParsingParams({TypeId::kInt, TypeId::kInt, TypeId::kMixed,
                           TypeId::kInt, TypeId::kInt, TypeId::kInt,
                           TypeId::kInt},
                          "neighbors10k.csv"),
        TypeParsingParams({TypeId::kDouble, TypeId::kDouble, TypeId::kDouble,
                           TypeId::kDouble, TypeId::kString},
                          "iris.csv", ',', false),
        TypeParsingParams({TypeId::kUndefined,TypeId::kUndefined,TypeId::kUndefined,
                           TypeId::kInt, TypeId::kString, TypeId::kDouble,
                           TypeId::kBigInt, TypeId::kMixed, TypeId::kBigInt,
                           TypeId::kMixed, TypeId::kInt},
                          "SimpleTypes.csv")));

TEST(TypeSystem, SumColumnDoubles) {
    std::vector<mo::TypedColumnData> col_data = CreateColumnData("iris.csv", ',', false);
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
