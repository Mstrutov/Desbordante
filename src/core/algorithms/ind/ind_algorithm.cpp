#include "ind_algorithm.h"

#include "config/names_and_descriptions.h"
#include "config/tabular_data/input_tables/option.h"

namespace algos {

INDAlgorithm::INDAlgorithm(std::vector<std::string_view> phase_names)
    : Algorithm(std::move(phase_names)) {
    RegisterOption(config::kTablesOpt(&input_tables_));
    MakeOptionsAvailable({config::kTablesOpt.GetName()});
}

void INDAlgorithm::LoadDataInternal() {
    schemas_ = std::make_shared<std::vector<RelationalSchema>>();
    for (auto const& input_table : input_tables_) {
        schemas_->emplace_back(input_table->GetRelationName(), input_table->GetColumnNames());
    }

    LoadINDAlgorithmDataInternal();
}

}  // namespace algos
