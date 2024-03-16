#include "algo_factory.h"

#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "algorithms/algorithms.h"
#include "algorithms/create_algorithm.h"
#include "algorithms/pipelines/typo_miner/typo_miner.h"
#include "config/names.h"
#include "tabular_data/input_tables_type.h"
#include "util/separator_validator.h"

namespace algos {

namespace {

boost::any GetAnyFromMap(StdParamsMap const& options, std::string_view option_name) {
    using std::is_same_v, std::decay;
    std::string const string_opt{option_name};
    auto it = options.find(string_opt);
    if (it == options.end()) {
        throw std::out_of_range("No option named \"" + string_opt + "\" in parameters.");
    }
    return it->second;
}

template <typename T>
T GetOptionValue(StdParamsMap const& options, std::string const& option_name) {
    return boost::any_cast<T>(GetAnyFromMap(options, option_name));
}

boost::any GetOrEmpty(StdParamsMap const& options, std::string_view option_name) {
    auto it = options.find(std::string{option_name});
    return it == options.end() ? boost::any{} : it->second;
}

}  // namespace

void ConfigureFromMap(Algorithm& algorithm, StdParamsMap const& options) {
    ConfigureFromFunction(algorithm, [&options](std::string_view option_name) {
        return GetOrEmpty(options, option_name);
    });
}

void LoadAlgorithm(Algorithm& algorithm, StdParamsMap const& options) {
    ConfigureFromFunction(algorithm, [&options](std::string_view option_name) {
        using namespace config::names;
        auto create_input_table = [](CSVConfig const& csv_config) -> config::InputTable {
            auto csv_parser = std::make_shared<CSVParser>(csv_config);
            csv_parser->ValidateSeparator();
            return csv_parser;
        };

        if (option_name == kTable && options.find(std::string{kTable}) == options.end()) {
            auto table = create_input_table(GetOptionValue<CSVConfig>(options, kCsvConfig));
            return boost::any{table};
        } else if (option_name == kTables && options.find(std::string{kTables}) == options.end()) {
            auto csv_configs = GetOptionValue<std::vector<CSVConfig>>(options, kCsvConfigs);
            if (csv_configs.empty()) {
                throw config::ConfigurationError("Expected collection of csv configs");
            }

            config::InputTables tables;
            tables.reserve(csv_configs.size());
            std::transform(csv_configs.begin(), csv_configs.end(), std::back_inserter(tables),
                           create_input_table);

            return boost::any{tables};
        }
        return GetOrEmpty(options, option_name);
    });
    algorithm.LoadData();
    ConfigureFromMap(algorithm, options);
}

std::unique_ptr<Algorithm> CreateAlgorithm(AlgorithmType algorithm_enum,
                                           StdParamsMap const& options) {
    std::unique_ptr<Algorithm> algorithm = CreateAlgorithmInstance(algorithm_enum);
    LoadAlgorithm(*algorithm, options);
    return algorithm;
}

std::unique_ptr<Algorithm> CreateTypoMiner(StdParamsMap const& options) {
    using config::names::kPreciseAlgorithm, config::names::kApproximateAlgorithm;
    AlgorithmType precise_algo = GetOptionValue<AlgorithmType>(options, kPreciseAlgorithm);
    AlgorithmType approx_algo = GetOptionValue<AlgorithmType>(options, kApproximateAlgorithm);
    std::unique_ptr<TypoMiner> typo_miner = std::make_unique<TypoMiner>(precise_algo, approx_algo);
    LoadAlgorithm(*typo_miner, options);
    return typo_miner;
}

std::unique_ptr<Algorithm> CreateAlgorithm(std::string const& algorithm_name,
                                           StdParamsMap const& options) {
    if (algorithm_name == "typo_miner") {
        return CreateTypoMiner(options);
    }
    AlgorithmType const algorithm_enum = AlgorithmType::_from_string_nocase(algorithm_name.c_str());
    return CreateAlgorithm(algorithm_enum, options);
}

}  // namespace algos
