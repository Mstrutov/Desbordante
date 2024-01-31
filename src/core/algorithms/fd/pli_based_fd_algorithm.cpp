#include "pli_based_fd_algorithm.h"

#include "config/equal_nulls/option.h"
#include "config/tabular_data/input_table/option.h"

namespace algos {

PliBasedFDAlgorithm::PliBasedFDAlgorithm(std::vector<std::string_view> phase_names)
    : FDAlgorithm(std::move(phase_names)) {
    RegisterOptions();
    MakeOptionsAvailable({config::TableOpt.GetName(), config::EqualNullsOpt.GetName()});
}

void PliBasedFDAlgorithm::RegisterOptions() {
    RegisterOption(config::TableOpt(&input_table_));
    RegisterOption(config::EqualNullsOpt(&is_null_equal_null_));
}

void PliBasedFDAlgorithm::LoadDataInternal() {
    relation_ = ColumnLayoutRelationData::CreateFrom(*input_table_, is_null_equal_null_);

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

void PliBasedFDAlgorithm::LoadData(std::shared_ptr<ColumnLayoutRelationData> data) {
    if (data->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty dataset: FD mining is meaningless.");
    }  // TODO: this has to be repeated for every "alternative" data load
    relation_ = std::move(data);
    ExecutePrepare();  // TODO: this has to be repeated for every "alternative" data load
}

}  // namespace algos
