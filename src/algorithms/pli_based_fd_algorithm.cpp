#include "algorithms/pli_based_fd_algorithm.h"

#include "util/config/names_and_descriptions.h"
#include "util/config/option_using.h"

namespace algos {

PliBasedFDAlgorithm::PliBasedFDAlgorithm(std::vector<std::string_view> phase_names,
                                         bool request_prepared_data)
    : FDAlgorithm(std::move(phase_names),
                  request_prepared_data ? util::config::ConfigurationStage::load_prepared_data
                                        : util::config::ConfigurationStage::load_data) {
    RegisterOptions();
}

void PliBasedFDAlgorithm::RegisterOptions() {
    DESBORDANTE_OPTION_USING;

    RegisterPrepLoadOption(Option{&relation_, kPliRelation, kDPliRelation});
}

void PliBasedFDAlgorithm::LoadDataInternal() {
    if (GetCurrentStage() == +util::config::ConfigurationStage::load_data) {
        assert(input_table_);
        relation_ = ColumnLayoutRelationData::CreateFrom(*input_table_, is_null_equal_null_);
    }

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
