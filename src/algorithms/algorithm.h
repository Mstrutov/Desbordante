#pragma once

#include <filesystem>
#include <string_view>
#include <typeindex>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/any.hpp>

#include "model/idataset_stream.h"
#include "parser/csv_parser.h"
#include "util/config/configuration.h"
#include "util/config/ioption.h"
#include "util/config/option.h"
#include "util/progress.h"

namespace algos {

class Algorithm {
private:
    util::Progress progress_;
    util::config::Configuration configuration_;

    // Clear the necessary fields for Execute to run repeatedly with different
    // configuration parameters on the same dataset.
    virtual void ResetState() = 0;

    virtual void LoadDataInternal() = 0;
    virtual unsigned long long ExecuteInternal() = 0;

protected:
    void AddProgress(double val) noexcept {
        progress_.AddProgress(val);
    }
    void SetProgress(double val) noexcept {
        progress_.SetProgress(val);
    }
    void ToNextProgressPhase() noexcept {
        progress_.ToNextProgressPhase();
    }

    util::config::ConfigurationStage GetCurrentStage() {
        return configuration_.GetCurrentStage();
    }

    void RegisterOption(util::config::IOption&& option) {
        configuration_.RegisterOption(std::move(option));
    }

    void RegisterPrepLoadOption(util::config::IOption&& option) {
        configuration_.RegisterOption(std::move(option),
                                      +util::config::ConfigurationStage::load_prepared_data);
    }

    void RegisterInitialLoadOption(util::config::IOption&& option) {
        configuration_.RegisterOption(std::move(option),
                                      +util::config::ConfigurationStage::load_data);
    }

    void RegisterInitialExecOption(util::config::IOption&& option) {
        configuration_.RegisterOption(std::move(option),
                                      +util::config::ConfigurationStage::execute);
    }

    // The constructor accepts vector of names of the mining algorithm phases.
    // NOTE: Pass an empty vector here if your algorithm does not have an implemented progress bar.
    explicit Algorithm(std::vector<std::string_view> phase_names,
                       util::config::ConfigurationStage starting_stage =
                               util::config::ConfigurationStage::load_data);

    explicit Algorithm(std::vector<std::string_view> phase_names,
                       util::config::Configuration::FuncTuple funcTuple);

public:
    constexpr static double kTotalProgressPercent = util::Progress::kTotalProgressPercent;

    Algorithm(Algorithm const& other) = delete;
    Algorithm& operator=(Algorithm const& other) = delete;
    Algorithm(Algorithm&& other) = delete;
    Algorithm& operator=(Algorithm&& other) = delete;
    virtual ~Algorithm() = default;

    void LoadData();

    unsigned long long Execute();

    void SetOption(std::string_view option_name, boost::any const& value = {}) {
        configuration_.SetOption(option_name, value);
    }

    std::pair<bool, std::string> SetOptionNoThrow(std::string_view option_name,
                                                  boost::any const& value = {}) {
        return configuration_.SetOptionNoThrow(option_name, value);
    }

    [[nodiscard]] std::unordered_set<std::string_view> GetNeededOptions() const {
        return configuration_.GetNeededOptions();
    }

    [[nodiscard]] std::unordered_set<std::string_view> GetPossibleOptions() const {
        return configuration_.GetPossibleOptions();
    }

    [[nodiscard]] std::string_view GetDescription(std::string_view option_name) const {
        return configuration_.GetDescription(option_name);
    }

    void UnsetOption(std::string_view option_name) noexcept {
        configuration_.UnsetOption(option_name);
    }

    [[nodiscard]] std::type_index GetTypeIndex(std::string_view option_name) const {
        return configuration_.GetTypeIndex(option_name);
    }

    [[nodiscard]] bool NeedsOption(std::string_view option_name) const {
        return configuration_.NeedsOption(option_name);
    }

    [[nodiscard]] bool OptionSettable(std::string_view option_name) const {
        return configuration_.OptionSettable(option_name);
    }

    [[nodiscard]] bool IsInitialAtStage(std::string_view option_name,
                                        util::config::ConfigurationStage stage) const {
        return configuration_.IsInitialAtStage(option_name, stage);
    }

    void ResetConfiguration() {
        configuration_.Reset();
    }

    // See util::Progress::GetProgress description
    std::pair<uint8_t, double> GetProgress() const noexcept {
        return progress_.GetProgress();
    }
    std::vector<std::string_view> const& GetPhaseNames() const noexcept {
        return progress_.GetPhaseNames();
    }
};

}  // namespace algos
