#pragma once

#include <array>
#include <functional>
#include <memory>
#include <string_view>
#include <typeindex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/any.hpp>

#include "util/config/configuration_stage.h"
#include "util/config/ioption.h"

namespace util::config {

class Configuration {
public:
    using OptionNameSet = std::unordered_set<std::string_view>;
    using AddExternalOptFunc = std::function<void(OptionNameSet&)>;
    using SetExternalOptFunc =
            std::function<std::pair<bool, std::string>(std::string_view, boost::any const&)>;
    using UnsetExternalOptFunc = std::function<void(std::string_view)>;
    using GetExternalTypeIndexFunc = std::function<std::type_index(std::string_view)>;
    using ResetExternalAlgoConfigFunc = std::function<void()>;
    using FuncTuple = std::tuple<AddExternalOptFunc, SetExternalOptFunc, UnsetExternalOptFunc,
                                 GetExternalTypeIndexFunc, ResetExternalAlgoConfigFunc>;

private:
    // All options the algorithm may use.
    std::unordered_map<std::string_view, std::unique_ptr<config::IOption>> possible_options_{};
    // All options that should be set for the algorithm to proceed (may change
    // as some options are set).
    OptionNameSet required_options_;
    // Maps an option's name to the names of options it made required.
    std::unordered_map<std::string_view, std::vector<std::string_view>> opt_parents_{};
    ConfigurationStage current_stage_ = ConfigurationStage::load_data;
    // Options that will become available after the corresponding stage is started.
    std::array<OptionNameSet, ConfigurationStage::_size()> initial_stage_options_;

    // Pipeline-related functions.
    // Adds names of options that are needed for algorithms inside a pipeline.
    AddExternalOptFunc add_external_needed_options_;
    // Tries to set value of an option of an algorithm inside a pipeline.
    SetExternalOptFunc set_external_option_;
    // Tries to unset value of an option of an algorithm inside a pipeline.
    UnsetExternalOptFunc unset_external_option_;
    // Gets the type index corresponding to the type of the option of an
    // algorithm inside a pipeline.
    GetExternalTypeIndexFunc get_external_type_index_;
    ResetExternalAlgoConfigFunc reset_external_algo_config_;

    void UnsetChildren(std::string_view parent_option_name) noexcept;

    void UnsetOptionInternal(std::string_view option_name) noexcept;

    // Add an option that becomes available as soon as a configuration stage
    // starts. This option is called "initial".
    void AddInitialOption(ConfigurationStage stage, std::string_view name);

    bool IsPipeline() const;

    Configuration(ConfigurationStage starting_stage, FuncTuple funcTuple);

public:
    explicit Configuration(ConfigurationStage starting_stage = ConfigurationStage::load_data);

    explicit Configuration(FuncTuple funcTuple);

    void RegisterOption(IOption&& option);

    void RegisterOption(IOption&& option, ConfigurationStage stage);

    // The pipeline should take care that its algorithms are called in the
    // correct way, so this method doesn't take its algorithms into account.
    void StartStage(ConfigurationStage stage);

    [[nodiscard]] ConfigurationStage GetCurrentStage() const;

    void Reset();

    void SetOption(std::string_view option_name, boost::any const& value = {});

    std::pair<bool, std::string> SetOptionNoThrow(std::string_view option_name,
                                                  boost::any const& value = {}) noexcept;

    void UnsetOption(std::string_view option_name) noexcept;

    // Get all options that are required and not set.
    [[nodiscard]] OptionNameSet GetNeededOptions() const;

    // Checks if the provided option has to be set to proceed. Does not account
    // for options from the algorithms in a pipeline.
    [[nodiscard]] bool NeedsOption(std::string_view option_name) const;

    [[nodiscard]] bool OptionSettable(std::string_view option_name) const;

    // Check if an option is initial (see above) at the provided configuration
    // stage. Does not account for the algorithms in a pipeline.
    [[nodiscard]] bool IsInitialAtStage(std::string_view option_name,
                                        ConfigurationStage stage) const;

    // In a normal algorithm: gets the type index of the option with the given name.
    // In a pipeline: gets the type index of the pipeline's option with the
    // given name if it is required, otherwise gets the type index its
    // algorithms need.
    [[nodiscard]] std::type_index GetTypeIndex(std::string_view option_name) const;

    [[nodiscard]] std::string_view GetDescription(std::string_view option_name) const;

    [[nodiscard]] std::unordered_set<std::string_view> GetPossibleOptions() const;
};

}  // namespace util::config
