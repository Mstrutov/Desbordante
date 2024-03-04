#include <gtest/gtest.h>

#include "algo_factory.h"
#include "all_csv_configs.h"
#include "config/names.h"
#include "fd/pfdtane/enums.h"
#include "fd/pfdtane/pfdtane.h"
#include "model/table/column_layout_relation_data.h"
#include "parser/csv_parser/csv_parser.h"

namespace tests {
namespace onam = config::names;

struct PFD {
    std::size_t lhs_id;
    std::size_t rhs_id;
    config::ErrorType expected_error;
};

struct PFDTaneMiningParams {
    algos::StdParamsMap params;
    unsigned int result_hash;

    PFDTaneMiningParams(unsigned int result_hash, config::ErrorType const error,
                        algos::ErrorMeasure error_measure, CSVConfig table_config)
        : params({{onam::kCsvConfig, table_config},
                  {onam::kError, error},
                  {onam::kErrorMeasure, error_measure}}),
          result_hash(result_hash) {}
};

struct PFDTaneValidationParams {
    std::vector<PFD> fds;
    algos::ErrorMeasure error_measure;
    CSVConfig table_config;
};

class TestPFDTaneMining : public ::testing::TestWithParam<PFDTaneMiningParams> {};

class TestPFDTaneValidation : public ::testing::TestWithParam<PFDTaneValidationParams> {};

TEST_P(TestPFDTaneMining, DefaultTest) {
    auto const &p = GetParam();
    auto mp = algos::StdParamsMap(p.params);
    auto algos = algos::CreateAndLoadAlgorithm<algos::PFDTane>(mp);
    algos->Execute();
    EXPECT_EQ(p.result_hash, algos->Fletcher16());
}

TEST_P(TestPFDTaneValidation, ErrorCalculationTest) {
    auto const &p = GetParam();
    double eps = 0.00001;
    auto table = std::make_shared<CSVParser>(p.table_config);
    auto relation = ColumnLayoutRelationData::CreateFrom(*table, true);
    for (auto const &[lhs_id, rhs_id, expected_error] : p.fds) {
        auto lhs = relation->GetColumnData(lhs_id).GetPositionListIndex();
        auto rhs = relation->GetColumnData(rhs_id).GetPositionListIndex();
        config::ErrorType error =
                algos::PFDTane::CalculateFdError(lhs, lhs->Intersect(rhs).get(), p.error_measure);
        EXPECT_NEAR(error, expected_error, eps);
    }
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(
        PFDTaneTestMiningSuite, TestPFDTaneMining,
        ::testing::Values(
            PFDTaneMiningParams(44381, 0.3, +algos::ErrorMeasure::per_value, kTestFD)
        ));

INSTANTIATE_TEST_SUITE_P(
        PFDTaneTestValidationSuite, TestPFDTaneValidation,
        ::testing::Values(
                    PFDTaneValidationParams({{2, 3, 0.0625}, {4, 5, 0.333333}, {3, 2, 0.291666}, {0, 1, 0.75},
                                             {1, 0, 0.0}, {4, 3, 0.099999}, {1, 5, 0.416666}, {5, 1, 0.0}}, +algos::ErrorMeasure::per_value, kTestFD),
                    PFDTaneValidationParams({{2, 3, 0.083333}, {4, 5, 0.333333}, {3, 2, 0.5}, {0, 1, 0.75},
                                             {1, 0, 0.0}, {4, 3, 0.083333}, {1, 5, 0.416666}, {5, 1, 0.0}}, +algos::ErrorMeasure::per_tuple, kTestFD)
                ));

// clang-format on

}  // namespace tests
