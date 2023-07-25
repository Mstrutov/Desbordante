#include "algorithms/algorithm.h"

#include <cassert>
#include <utility>

namespace algos {

Algorithm::Algorithm(std::vector<std::string_view> phase_names,
                     util::config::ConfigurationStage starting_stage)
    : progress_(std::move(phase_names)), configuration_(starting_stage) {}

Algorithm::Algorithm(std::vector<std::string_view> phase_names,
                     util::config::Configuration::FuncTuple funcTuple)
    : progress_(std::move(phase_names)), configuration_(std::move(funcTuple)) {}

void Algorithm::LoadData() {
    if (configuration_.GetCurrentStage() == +util::config::ConfigurationStage::execute) {
        throw std::logic_error("Data has already been processed.");
    }
    if (!GetNeededOptions().empty()) {
        throw std::logic_error("All options need to be set before starting processing.");
    }
    LoadDataInternal();
    configuration_.StartStage(util::config::ConfigurationStage::execute);
}

unsigned long long Algorithm::Execute() {
    if (configuration_.GetCurrentStage() != +util::config::ConfigurationStage::execute) {
        throw std::logic_error("Data must be processed before execution.");
    }
    if (!GetNeededOptions().empty()) {
        throw std::logic_error("All options need to be set before execution.");
    }
    progress_.ResetProgress();
    ResetState();
    auto time_ms = ExecuteInternal();
    ResetConfiguration();
    return time_ms;
}

}  // namespace algos

