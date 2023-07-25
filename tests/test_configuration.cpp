#include <gtest/gtest.h>

#include "util/config/configuration.h"
#include "util/config/option.h"

namespace tests {

using util::config::Configuration;
using util::config::ConfigurationStage;
using util::config::Option;

using FuncTuple = Configuration::FuncTuple;
using OptionNameSet = Configuration::OptionNameSet;
using AddExternalOptFunc = Configuration::AddExternalOptFunc;
using SetExternalOptFunc = Configuration::SetExternalOptFunc;
using UnsetExternalOptFunc = Configuration::UnsetExternalOptFunc;
using GetExternalTypeIndexFunc = Configuration::GetExternalTypeIndexFunc;
using ResetExternalAlgoConfigFunc = Configuration::ResetExternalAlgoConfigFunc;

class ConfigurationTesting : public ::testing::Test {
protected:
    Configuration config;
    ConfigurationStage load_data = ConfigurationStage::load_data;
    ConfigurationStage execute = ConfigurationStage::execute;
    ConfigurationStage load_prepared_data = ConfigurationStage::load_prepared_data;

    void RegisterInitialLoadOption(util::config::IOption&& option) {
        config.RegisterOption(std::move(option), load_data);
    }

    void RegisterInitialExecOption(util::config::IOption&& option) {
        config.RegisterOption(std::move(option), execute);
    }

    void RegisterOption(util::config::IOption&& option) {
        config.RegisterOption(std::move(option));
    }
};

TEST_F(ConfigurationTesting, TestEmpty) {
    ASSERT_EQ(config.GetCurrentStage(), load_data);
    ASSERT_NO_THROW(config.UnsetOption("whatever"));
    ASSERT_EQ(config.GetNeededOptions(), OptionNameSet{});
    ASSERT_FALSE(config.NeedsOption("whatever"));
    ASSERT_FALSE(config.IsInitialAtStage("whatever", execute));
    ASSERT_EQ(config.GetTypeIndex("whatever"), typeid(void));
    ASSERT_NO_THROW(config.Reset());
    ASSERT_NO_THROW(config.StartStage(execute));
    ASSERT_EQ(config.GetCurrentStage(), execute);
    ASSERT_NO_THROW(config.Reset());
    ASSERT_EQ(config.GetNeededOptions(), OptionNameSet{});
    ASSERT_THROW(config.SetOption("whatever"), std::invalid_argument);
    ASSERT_THROW(config.SetOption("whatever", 2), std::invalid_argument);
}

TEST_F(ConfigurationTesting, BasicConfiguration) {
    int var1;
    const int var1val = 4;
    double var2;
    const double var2default = 4.1;
    char var3;
    constexpr char var3val1 = ',';
    constexpr char var3val2 = ';';

    config.RegisterOption(Option{&var1, "name1", "desc1"}, load_data);
    config.RegisterOption(Option{&var2, "name2", "desc2", var2default}, load_data);
    config.RegisterOption(Option{&var3, "name3", "desc3", var3val1}, execute);

    ASSERT_EQ(config.GetNeededOptions(), (OptionNameSet{"name1", "name2"}));

    ASSERT_THROW(config.SetOption("name1"), std::invalid_argument);
    ASSERT_NO_THROW(config.SetOption("name1", var1val));
    ASSERT_EQ(var1, var1val);

    ASSERT_EQ(config.GetNeededOptions(), OptionNameSet{"name2"});

    ASSERT_THROW(config.SetOption("name2", 1), std::invalid_argument);
    ASSERT_NO_THROW(config.SetOption("name2"));
    ASSERT_EQ(var2, var2default);

    ASSERT_EQ(config.GetNeededOptions(), OptionNameSet{});

    config.Reset();

    ASSERT_EQ(config.GetNeededOptions(), (OptionNameSet{"name1", "name2"}));
    ASSERT_NO_THROW(config.SetOption("name1", var1val));
    ASSERT_NO_THROW(config.SetOption("name2"));

    ASSERT_EQ(config.GetNeededOptions(), OptionNameSet{});

    ASSERT_THROW(config.SetOption("name3"), std::invalid_argument);

    config.StartStage(execute);
    ASSERT_EQ(config.GetNeededOptions(), OptionNameSet{"name3"});

    ASSERT_NO_THROW(config.SetOption("name3"));
    ASSERT_EQ(var3, var3val1);

    ASSERT_NO_THROW(config.SetOption("name3", var3val2));
    ASSERT_EQ(var3, var3val2);

    ASSERT_EQ(config.GetNeededOptions(), OptionNameSet{});

    ASSERT_NO_THROW(config.UnsetOption("name3"));
    ASSERT_EQ(config.GetNeededOptions(), OptionNameSet{"name3"});
}

TEST_F(ConfigurationTesting, ConditionalConfiguration) {
    using IndicesType = std::vector<unsigned int>;
    enum class Metric {
        cosine,
        euclidean,
        levenshtein,
    };
    enum class MetricAlgo {
        brute,
        calipers,
        approx,
    };

    bool eq_nulls;
    Metric metric;
    IndicesType rhs_indices;
    MetricAlgo algo;
    unsigned int q;
    double parameter;

    auto check_parameter = [](long double parameter) {
        if (parameter < 0) throw std::invalid_argument("Parameter out of range");
    };
    auto ind_norm = [](IndicesType& indices) {
        std::sort(indices.begin(), indices.end());
        indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
    };
    auto check_rhs = [](IndicesType const& rhs_indices) {
        if (rhs_indices.back() >= 3) {
            throw std::invalid_argument(
                    "Column index should be less than the number of columns in the dataset.");
        }
    };
    auto need_algo_and_q = [&metric]([[maybe_unused]] IndicesType const& _) {
        return metric == Metric::cosine;
    };
    auto need_algo_only = [&metric](IndicesType const& rhs_indices) {
        return metric == Metric::levenshtein || rhs_indices.size() != 1;
    };
    auto algo_check = [&metric, &rhs_indices](MetricAlgo metric_algo) {
        if (metric_algo == MetricAlgo::calipers) {
            if (!(metric == Metric::euclidean && rhs_indices.size() == 2))
                throw std::invalid_argument(
                        "\"calipers\" algorithm is only available for 2-dimensional RHS and "
                        "\"euclidean\" metric.");
        }
    };
    auto q_check = [](unsigned int q) {
        if (q <= 0) throw std::invalid_argument("Q-gram length should be greater than zero.");
    };

    RegisterInitialLoadOption(Option{&eq_nulls, "kEqualNulls", "kDEqualNulls", true});
    RegisterInitialExecOption(
            Option{&parameter, "kParameter", "kDParameter"}.SetValueCheck(check_parameter));
    RegisterInitialExecOption(
            Option{&metric, "kMetric", "kDMetric"}.SetConditionalOpts({{{}, {"kRhsIndices"}}}));
    RegisterOption(
            Option{&algo, "kMetricAlgorithm", "kDMetricAlgorithm"}.SetValueCheck(algo_check));
    RegisterOption(Option{&q, "kQGramLength", "kDQGramLength", 2u}.SetValueCheck(q_check));
    RegisterOption(
            Option{&rhs_indices, "kRhsIndices", "kDRhsIndices"}
                    .SetNormalizeFunc(ind_norm)
                    .SetValueCheck(check_rhs)
                    .SetConditionalOpts({{need_algo_and_q, {"kMetricAlgorithm", "kQGramLength"}},
                                         {need_algo_only, {"kMetricAlgorithm"}}}));

    ASSERT_EQ(config.GetTypeIndex("kRhsIndices"), typeid(IndicesType));

    ASSERT_FALSE(config.GetNeededOptions().empty());
    ASSERT_NO_THROW(config.SetOption("kEqualNulls"));
    ASSERT_TRUE(config.GetNeededOptions().empty());

    config.StartStage(execute);

    ASSERT_EQ(config.GetNeededOptions(), (OptionNameSet{"kParameter", "kMetric"}));
    ASSERT_THROW(config.SetOption("kParameter"), std::invalid_argument);
    ASSERT_THROW(config.SetOption("kParameter", -2.0), std::invalid_argument);
    ASSERT_NO_THROW(config.SetOption("kParameter", 0.2));

    ASSERT_NO_THROW(config.SetOption("kMetric", Metric::cosine));
    ASSERT_THROW(config.SetOption("kRhsIndices", IndicesType{3, 0}), std::invalid_argument);
    ASSERT_NO_THROW(config.SetOption("kRhsIndices", IndicesType{0, 2}));
    OptionNameSet needed_cosine = config.GetNeededOptions();

    ASSERT_TRUE(needed_cosine.find("kQGramLength") != needed_cosine.end());
    ASSERT_EQ(needed_cosine, (OptionNameSet{"kMetricAlgorithm", "kQGramLength"}));
    ASSERT_THROW(config.SetOption("kQGramLength", 0u), std::invalid_argument);

    ASSERT_NO_THROW(config.SetOption("kMetric", Metric::levenshtein));
    ASSERT_EQ(config.GetNeededOptions(), (OptionNameSet{"kRhsIndices"}));
    ASSERT_NO_THROW(config.SetOption("kRhsIndices", IndicesType{0}));
    OptionNameSet needed_levenshtein = config.GetNeededOptions();
    ASSERT_TRUE(needed_levenshtein.find("kQGramLength") == needed_levenshtein.end());
    ASSERT_EQ(needed_levenshtein, (OptionNameSet{"kMetricAlgorithm"}));
    ASSERT_NO_THROW(config.SetOption("kMetricAlgorithm", MetricAlgo::brute));
    ASSERT_TRUE(config.GetNeededOptions().empty());

    config.Reset();
    ASSERT_EQ(config.GetNeededOptions(), (OptionNameSet{"kParameter", "kMetric"}));
    ASSERT_NO_THROW(config.SetOption("kParameter", 0.2));
    ASSERT_NO_THROW(config.SetOption("kMetric", Metric::euclidean));
    ASSERT_EQ(config.GetNeededOptions(), (OptionNameSet{"kRhsIndices"}));
    ASSERT_NO_THROW(config.SetOption("kRhsIndices", IndicesType{0}));
    ASSERT_TRUE(config.GetNeededOptions().empty());
    ASSERT_THROW(config.SetOption("kMetricAlgorithm", MetricAlgo::calipers), std::invalid_argument);
    ASSERT_NO_THROW(config.SetOption("kRhsIndices", IndicesType{0, 1, 2}));
    ASSERT_EQ(config.GetNeededOptions(), (OptionNameSet{"kMetricAlgorithm"}));
    ASSERT_THROW(config.SetOption("kMetricAlgorithm", MetricAlgo::calipers), std::invalid_argument);
    ASSERT_NO_THROW(config.SetOption("kRhsIndices", IndicesType{0, 1}));
    ASSERT_EQ(config.GetNeededOptions(), (OptionNameSet{"kMetricAlgorithm"}));
    ASSERT_NO_THROW(config.SetOption("kMetricAlgorithm", MetricAlgo::calipers));
    ASSERT_TRUE(config.GetNeededOptions().empty());
}

TEST_F(ConfigurationTesting, PipelineConfiguration) {
    static std::type_index const void_index = typeid(void);

    Configuration config1;
    Configuration config2;
    std::unique_ptr<Configuration> pipeline_configuration;

    auto try_set_option = [&config1, &config2](
                                  std::string_view option_name, boost::any const& value1,
                                  boost::any const& value2) -> std::pair<bool, std::string> {
        bool succeeded = false;
        std::string error_text;
        for (auto [algo, value] :
             {std::make_pair(&config1, value1), std::make_pair(&config2, value2)}) {
            if (!algo->OptionSettable(option_name)) continue;
            auto [success_config, error_text_config] = algo->SetOptionNoThrow(option_name, value);
            if (success_config) succeeded = true;
            if (succeeded) continue;
            error_text = error_text_config;
        }
        if (succeeded) error_text = "";
        return {succeeded, error_text};
    };

    AddExternalOptFunc add_external_needed_opts = [&config1,
                                                   &config2](OptionNameSet& needed_options) {
        OptionNameSet opts1 = config1.GetNeededOptions();
        OptionNameSet opts2 = config2.GetNeededOptions();
        needed_options.insert(opts1.begin(), opts1.end());
        needed_options.insert(opts2.begin(), opts2.end());
    };
    SetExternalOptFunc set_external_opt = [try_set_option](std::string_view option_name,
                                                           boost::any const& value) {
        if (option_name == "kError") {
            if (value.type() != typeid(char)) {
                throw std::invalid_argument("Incorrect type.");
            }
            if (value.empty()) {
                throw std::invalid_argument("Must specify value.");
            }
            auto error = boost::any_cast<char>(value);
            if (error == 0) {
                throw std::invalid_argument("0 is meaningless");
            }
            return try_set_option(option_name, char{0}, value);
        }
        return try_set_option(option_name, value, value);
    };
    UnsetExternalOptFunc unset_external_opt = [&config1, &config2](std::string_view option_name) {
        config1.UnsetOption(option_name);
        config2.UnsetOption(option_name);
    };
    GetExternalTypeIndexFunc get_external_type_index = [&config1,
                                                        &config2](std::string_view option_name) {
        std::type_index idx1 = config1.GetTypeIndex(option_name);
        std::type_index idx2 = config2.GetTypeIndex(option_name);
        if (idx1 != void_index) {
            if (idx2 != void_index && !config1.NeedsOption(option_name)) return idx2;
            return idx1;
        }
        return idx2;
    };
    ResetExternalAlgoConfigFunc reset_external_config = [&config1, &config2]() {
        config1.Reset();
        config2.Reset();
    };

    pipeline_configuration = std::make_unique<Configuration>(
            FuncTuple{add_external_needed_opts, set_external_opt, unset_external_opt,
                      get_external_type_index, reset_external_config});

    int opt11;
    int opt12;
    std::string opt21;
    std::string opt22;
    std::string opt2p;
    char opt3;
    unsigned long long opt4;

    auto opt3_check = [](char opt_val) {
        if (opt_val == '\0') {
            throw std::invalid_argument("Can't be null!");
        }
    };

    config1.RegisterOption(Option{&opt11, "opt1", "dopt1"}, ConfigurationStage::load_data);
    config1.RegisterOption(Option{&opt21, "opt2", "dopt2"}, ConfigurationStage::load_data);
    config2.RegisterOption(Option{&opt12, "opt1", "dopt1"}, ConfigurationStage::load_data);
    config2.RegisterOption(Option{&opt22, "opt2", "dopt2"}, ConfigurationStage::load_data);
    pipeline_configuration->RegisterOption(Option{&opt2p, "opt2", "dopt2"},
                                           ConfigurationStage::load_data);
    config1.RegisterOption(Option{&opt3, "opt3", "dopt3", '2'}.SetValueCheck(opt3_check),
                           ConfigurationStage::load_data);
    config2.RegisterOption(Option{&opt4, "opt4", "dopt4"}, ConfigurationStage::load_data);

    ASSERT_EQ(pipeline_configuration->GetNeededOptions(),
              (OptionNameSet{"opt1", "opt2", "opt3", "opt4"}));
    ASSERT_NO_THROW(pipeline_configuration->SetOption("opt1", 1));
    ASSERT_EQ(opt11, 1);
    ASSERT_EQ(opt12, 1);
    ASSERT_EQ(pipeline_configuration->GetNeededOptions(), (OptionNameSet{"opt2", "opt3", "opt4"}));
    ASSERT_NO_THROW(pipeline_configuration->UnsetOption("opt1"));
    ASSERT_EQ(pipeline_configuration->GetNeededOptions(),
              (OptionNameSet{"opt1", "opt2", "opt3", "opt4"}));

    ASSERT_NO_THROW(pipeline_configuration->SetOption("opt1", 1));
    ASSERT_EQ(opt11, 1);
    ASSERT_EQ(opt12, 1);
    ASSERT_NO_THROW(pipeline_configuration->SetOption("opt2", std::string("text")));
    ASSERT_EQ(opt21, std::string("text"));
    ASSERT_EQ(opt22, std::string("text"));
    ASSERT_EQ(opt2p, std::string("text"));
    ASSERT_EQ(pipeline_configuration->GetNeededOptions(), (OptionNameSet{"opt3", "opt4"}));
    ASSERT_THROW(pipeline_configuration->SetOption("opt3", '\0'), std::invalid_argument);
    ASSERT_NO_THROW(pipeline_configuration->SetOption("opt3"));
    ASSERT_EQ(opt3, '2');
    ASSERT_NO_THROW(pipeline_configuration->SetOption("opt3", ';'));
    ASSERT_EQ(opt3, ';');
    ASSERT_THROW(pipeline_configuration->SetOption("opt5", 11), std::invalid_argument);
    ASSERT_THROW(pipeline_configuration->SetOption("opt4", 7), std::invalid_argument);
    ASSERT_NO_THROW(pipeline_configuration->SetOption("opt4", 3ull));
    ASSERT_EQ(opt4, 3);
    ASSERT_TRUE(pipeline_configuration->GetNeededOptions().empty());
    ASSERT_NO_THROW(pipeline_configuration->Reset());
    ASSERT_EQ(pipeline_configuration->GetNeededOptions(),
              (OptionNameSet{"opt1", "opt2", "opt3", "opt4"}));

    ASSERT_EQ(pipeline_configuration->GetTypeIndex("opt3"), std::type_index(typeid(char)));
    ASSERT_EQ(pipeline_configuration->GetTypeIndex("opt2"), std::type_index(typeid(std::string)));
}

}  // namespace tests
