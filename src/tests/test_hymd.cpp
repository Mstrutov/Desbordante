#include <cstddef>
#include <filesystem>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/preprocessing/similarity_measure/levenshtein_similarity_measure.h"
#include "algorithms/md/hymd/utility/md_less.h"
#include "all_csv_configs.h"
#include "config/names.h"
#include "config/tabular_data/input_table_type.h"
#include "model/index.h"
#include "parser/csv_parser/csv_parser.h"

namespace {
auto GetCardinality(std::vector<model::md::DecisionBoundary> const& lhs_bounds) {
    return std::count_if(lhs_bounds.begin(), lhs_bounds.end(),
                         [](auto bound) { return bound != 0.0; });
}
}  // namespace

namespace tests {

class HyMDTest : public ::testing::Test {
protected:
    static algos::StdParamsMap GetParamMap(CSVConfig const& csv_config,
                                           std::optional<std::size_t> min_support = std::nullopt,
                                           bool prune_nondisjoint = true,
                                           model::md::DecisionBoundary minimum_similarity = 0.7) {
        using namespace config::names;
        using namespace algos::hymd;
        config::InputTable table = std::make_unique<CSVParser>(csv_config);
        std::shared_ptr<SimilarityMeasureCreator> measure_creator = std::make_shared<
                preprocessing::similarity_measure::LevenshteinSimilarityMeasure::Creator>(
                minimum_similarity);
        std::vector<std::tuple<std::string, std::string, std::shared_ptr<SimilarityMeasureCreator>>>
                column_matches_option;
        std::size_t const number_of_columns = table->GetNumberOfColumns();
        column_matches_option.reserve(number_of_columns);
        for (model::Index i = 0; i < number_of_columns; ++i) {
            std::string const column_name = table->GetColumnName(i);
            column_matches_option.emplace_back(column_name, column_name, measure_creator);
        }
        algos::StdParamsMap param_map = {
                {kLeftTable, table},
                {kPruneNonDisjoint, prune_nondisjoint},
                {kColumnMatches, column_matches_option},
        };
        if (min_support) {
            param_map[kMinSupport] = *min_support;
        }

        return param_map;
    }
};

TEST_F(HyMDTest, DrunkAnimalsNormal) {
    using model::md::DecisionBoundary, model::Index, model::MD;
    using MdPair = std::pair<std::vector<DecisionBoundary>, std::pair<Index, DecisionBoundary>>;
    auto param_map = GetParamMap(kDrunkAnimals);
    auto hymd = algos::CreateAndLoadAlgorithm<algos::hymd::HyMD>(param_map);
    std::vector<MdPair> expected = {
            {{0.0, 0.0, 0.0, 0.75}, {2, 0.75}},
            {{0.0, 0.0, 0.75, 0.0}, {3, 0.75}},
    };
    std::vector<MdPair> actual;
    algos::ConfigureFromMap(*hymd, param_map);
    hymd->Execute();
    auto actual_mds = hymd->MdList();
    std::transform(actual_mds.begin(), actual_mds.end(), std::back_inserter(actual),
                   [](model::MD const& md) -> MdPair {
                       return {md.GetLhsDecisionBounds(), md.GetRhs()};
                   });
    ASSERT_EQ(expected, actual);
}

TEST_F(HyMDTest, DrunkAnimalsNoLimits) {
    using model::md::DecisionBoundary, model::Index, model::MD, algos::hymd::utility::MdPair;
    auto param_map = GetParamMap(kDrunkAnimals, 0, false, 0.0);
    auto hymd = algos::CreateAndLoadAlgorithm<algos::hymd::HyMD>(param_map);
    std::vector<MdPair> expected = {
            {{0, 0, 0, 0}, {1, 1 / 6.}},   {{0.2, 0, 0, 0}, {3, 1}},
            {{0.2, 0, 0, 0}, {2, 0.75}},   {{1, 0, 0, 0}, {1, 1}},
            {{1, 0, 0, 0}, {2, 1}},        {{0, 0, 0.75, 0}, {3, 0.75}},
            {{0, 0, 0, 0.75}, {2, 0.75}},  {{0, 0, 0, 1}, {0, 1 / 5.}},

            {{0, 0, 1, 1}, {0, 1}},        {{0, 0, 1, 1}, {1, 1}},
            {{0, 1, 0, 0.75}, {0, 1}},     {{0, 1, 0, 0.75}, {2, 1}},
            {{0, 1, 0, 0.75}, {3, 1}},     {{0, 1, 0.75, 0}, {0, 1}},
            {{0, 1, 0.75, 0}, {2, 1}},     {{0, 1, 0.75, 0}, {3, 1}},
            {{0.125, 0, 0, 0.75}, {3, 1}}, {{0.125, 0, 0, 0.75}, {0, 1 / 5.}},
            {{0.125, 0, 0.75, 0}, {3, 1}}, {{0.125, 0, 0.75, 0}, {0, 1 / 5.}},
            {{0.125, 0, 1, 0}, {0, 1}},    {{0.125, 0, 1, 0}, {1, 1}},
            {{0.125, 1, 0, 0}, {0, 1}},    {{0.125, 1, 0, 0}, {2, 1}},
            {{0.125, 1, 0, 0}, {3, 1}},
    };
    std::sort(expected.begin(), expected.end(), algos::hymd::utility::MdLessPairs);
    std::vector<MdPair> actual;
    algos::ConfigureFromMap(*hymd, param_map);
    hymd->Execute();
    auto actual_mds = hymd->MdList();
    std::transform(actual_mds.begin(), actual_mds.end(), std::back_inserter(actual),
                   [](model::MD const& md) -> MdPair {
                       return {md.GetLhsDecisionBounds(), md.GetRhs()};
                   });
    ASSERT_EQ(expected, actual);
}

TEST_F(HyMDTest, AdultMDs) {
    using model::md::DecisionBoundary, model::Index, model::MD;
    using IndexPair = std::pair<Index, Index>;
    auto param_map = GetParamMap(kAdult);
    auto hymd = algos::CreateAndLoadAlgorithm<algos::hymd::HyMD>(param_map);
    algos::ConfigureFromMap(*hymd, param_map);
    hymd->Execute();
    auto actual_mds = hymd->MdList();
    for (auto [left, right] : {IndexPair{4, 3}, IndexPair{3, 4}}) {
        auto it = std::find_if(actual_mds.begin(), actual_mds.end(), [left, right](MD const& md) {
            Index rhs_index = md.GetRhs().first;
            if (rhs_index != right) return false;
            auto lhs_bounds = md.GetLhsDecisionBounds();
            auto card = GetCardinality(lhs_bounds);
            if (card != 1) return false;
            return lhs_bounds.at(left) == 1.0;
        });
        if (it == actual_mds.end()) {
            std::stringstream found_mds_string;
            for (MD const& md : actual_mds) {
                found_mds_string << md.ToStringShort() << ",";
            }
            found_mds_string.seekp(-1, std::stringstream::cur);
            FAIL() << "No FD-like MD found (left = " << left << ", right = " << right
                   << "). Found MDs: " << found_mds_string.str();
        }
    }
    ASSERT_EQ(111u, actual_mds.size());
}

}  // namespace tests
