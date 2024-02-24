#include "pli_based_fd_algorithm.h"

#include "config/equal_nulls/option.h"
#include "config/tabular_data/input_table/option.h"

namespace algos {

PliBasedFDAlgorithm::PliBasedFDAlgorithm(
        std::vector<std::string_view> phase_names,
        std::optional<ColumnLayoutRelationDataManager> relation_manager)
    : FDAlgorithm(std::move(phase_names)),
      relation_manager_(relation_manager.has_value()
                                ? *relation_manager
                                : ColumnLayoutRelationDataManager{
                                          &input_table_, &is_null_equal_null_, &relation_}) {
    if (relation_manager.has_value()) return;
    RegisterRelationManagerOptions();
    MakeOptionsAvailable({config::TableOpt.GetName(), config::EqualNullsOpt.GetName()});
}

void PliBasedFDAlgorithm::RegisterRelationManagerOptions() {
    RegisterOption(config::TableOpt(&input_table_));
    RegisterOption(config::EqualNullsOpt(&is_null_equal_null_));
}

void PliBasedFDAlgorithm::LoadDataInternal() {
    relation_ = relation_manager_.GetRelation();

    if (relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty dataset: FD mining is meaningless.");
    }
}

std::vector<Column const*> PliBasedFDAlgorithm::GetKeys() const {
    assert(relation_ != nullptr);

    std::vector<Column const*> keys;
    for (ColumnData const& col : relation_->GetColumnData()) {
        if (col.GetPositionListIndex()->AllValuesAreUnique()) {
            keys.push_back(col.GetColumn());
        }
    }

    return keys;
}

}  // namespace algos
