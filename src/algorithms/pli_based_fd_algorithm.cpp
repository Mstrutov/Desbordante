#include "algorithms/pli_based_fd_algorithm.h"

namespace algos {

PliBasedFDAlgorithm::PliBasedFDAlgorithm(std::vector<std::string_view> phase_names)
        : FDAlgorithm(std::move(phase_names)) {}

void PliBasedFDAlgorithm::FitFd(model::IDatasetStream& data_stream) {
    relation_ = ColumnLayoutRelationData::CreateFrom(data_stream, is_null_equal_null_);

    if (relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty dataset: FD mining is meaningless.");
    }
}

std::vector<Column const*> PliBasedFDAlgorithm::GetKeys() const {
    assert(relation_ != nullptr);

    std::vector<Column const*> keys;
    for (ColumnData const& col : relation_->GetColumnData()) {
        if (col.GetPositionListIndex()->GetNumNonSingletonCluster() == 0) {
            keys.push_back(col.GetColumn());
        }
    }

    return keys;
}

void PliBasedFDAlgorithm::Fit(std::shared_ptr<ColumnLayoutRelationData> data) {
    // TODO: this has to be repeated for every "alternative" Fit
    if (configuration_.GetCurrentStage() != +config::ConfigurationStage::fit)
        throw std::logic_error("Incorrect algorithm execution order: Fit.");
    if (data->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty dataset: FD mining is meaningless.");
    }
    number_of_columns_ = data->GetNumColumns();
    relation_ = std::move(data);
    // TODO: this has to be repeated for every "alternative" Fit
    configuration_.StartStage(config::ConfigurationStage::execute);
}

}  // namespace algos
