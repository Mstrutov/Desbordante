#pragma once

#include <type_traits>

#include <boost/any.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/program_options.hpp>
#include <enum.h>

#include "algorithms/algorithms.h"
#include "algorithms/create_primitive.h"
#include "algorithms/legacy_algorithms.h"
#include "algorithms/options/names.h"
#include "algorithms/typo_miner.h"
#include "algorithms/cfd/cfd_discovery.h"

namespace algos {

using StdParamsMap = std::unordered_map<std::string, boost::any>;

namespace details {

template <typename OptionMap>
boost::any ExtractAnyFromMap(OptionMap&& options, std::string_view const& option_name) {
    using std::is_same_v, std::decay, boost::program_options::variables_map;
    const std::string string_opt{option_name};
    auto it = options.find(string_opt);
    if (it == options.end()) {
        throw std::out_of_range("No option named \"" + string_opt + "\" in parameters.");
    }
    if constexpr (is_same_v<typename decay<OptionMap>::type, variables_map>) {
        return options.extract(it).mapped().value();
    } else {
        return options.extract(it).mapped();
    }
}

template <typename T, typename OptionMap>
T ExtractOptionValue(OptionMap&& options, std::string const& option_name) {
    return boost::any_cast<T>(ExtractAnyFromMap(std::forward<OptionMap>(options), option_name));
}

template <typename ParamsMap>
ACAlgorithm::Config CreateAcAlgorithmConfigFromMap(ParamsMap params) {
    namespace onam = config::names;
    ACAlgorithm::Config c;

    c.data = std::filesystem::current_path() / "input_data" /
             ExtractOptionValue<std::string>(params, onam::kData);
    c.separator = ExtractOptionValue<char>(params, onam::kSeparator);
    c.has_header = ExtractOptionValue<bool>(params, onam::kHasHeader);
    c.bin_operation = ExtractOptionValue<char>(params, onam::kBinaryOperation);
    c.fuzziness = ExtractOptionValue<double>(params, onam::kFuzziness);
    if (c.fuzziness <= 0 || c.fuzziness > 1) {
        throw std::invalid_argument(
                "Fuzziness value must belong to the interval: (0, 1]");
    }
    c.p_fuzz = ExtractOptionValue<double>(params, onam::kFuzzinessProbability);
    if (c.p_fuzz <= 0 || c.p_fuzz > 1) {
        throw std::invalid_argument(
                "FuzzinessProbability value must belong to the interval: (0, 1]");
    }
    c.weight = ExtractOptionValue<double>(params, onam::kWeight);
    if (c.weight <= 0 || c.weight > 1) {
        throw std::invalid_argument("Weight value must belong to the interval: (0, 1]");
    }
    c.bumps_limit = ExtractOptionValue<size_t>(params, onam::kBumpsLimit);
    c.iterations_limit = ExtractOptionValue<size_t>(params, onam::kIterationsLimit);
    if (c.iterations_limit < 1) {
        throw std::invalid_argument("IterationsLimit value should not be less than one");
    }
    c.pairing_rule = ExtractOptionValue<std::string>(params, onam::kPairingRule);

    return c;
}

template <typename ParamsMap>
std::unique_ptr<Primitive> CreateAcAlgorithmInstance(ParamsMap&& params) {
    ACAlgorithm::Config const config =
        CreateAcAlgorithmConfigFromMap(std::forward<ParamsMap>(params));
    return std::make_unique<ACAlgorithm>(config);
}

}  // namespace details

template <typename OptionMap>
void ConfigureFromMap(Primitive& primitive, OptionMap&& options) {
    std::unordered_set<std::string_view> needed;
    while (!(needed = primitive.GetNeededOptions()).empty()) {
        for (std::string_view const& option_name : needed) {
            if (options.find(std::string{option_name}) == options.end()) {
                primitive.SetOption(option_name);
            } else {
                primitive.SetOption(option_name, details::ExtractAnyFromMap(options, option_name));
            }
        }
    }
}

template <typename OptionMap>
void LoadPrimitive(Primitive& prim, OptionMap&& options) {
    ConfigureFromMap(prim, options);
    auto parser = CSVParser{
            details::ExtractOptionValue<std::string>(options, config::names::kData),
            details::ExtractOptionValue<char>(options, config::names::kSeparator),
            details::ExtractOptionValue<bool>(options, config::names::kHasHeader)};
    prim.Fit(parser);
    ConfigureFromMap(prim, options);
}

template <typename T, typename OptionMap>
std::unique_ptr<T> CreateAndLoadPrimitive(OptionMap&& options) {
    std::unique_ptr<T> prim = std::make_unique<T>();
    LoadPrimitive(*prim, std::forward<OptionMap>(options));
    return prim;
}

template <typename OptionMap>
std::unique_ptr<Primitive> CreatePrimitive(PrimitiveType primitive_enum, OptionMap&& options) {
    std::unique_ptr<Primitive> primitive = CreatePrimitiveInstance(primitive_enum);
    LoadPrimitive(*primitive, std::forward<OptionMap>(options));
    return primitive;
}

template <typename OptionMap>
std::unique_ptr<Primitive> CreateTypoMiner(OptionMap&& options) {
    PrimitiveType precise_algo = details::ExtractOptionValue<PrimitiveType>(
            options, config::names::kPreciseAlgorithm);
    PrimitiveType approx_algo = details::ExtractOptionValue<PrimitiveType>(
            options, config::names::kApproximateAlgorithm);
    std::unique_ptr<TypoMiner> typo_miner = std::make_unique<TypoMiner>(precise_algo, approx_algo);
    LoadPrimitive(*typo_miner, std::forward<OptionMap>(options));
    return typo_miner;
}

template <typename OptionMap>
std::unique_ptr<Primitive> CreatePrimitive(std::string const& primitive_name,
                                           OptionMap&& options) {
    if (primitive_name == "ac") {
        return details::CreateAcAlgorithmInstance(std::forward<OptionMap>(options));
    }
    if (primitive_name == "typo_miner") {
        return CreateTypoMiner(std::forward<OptionMap>(options));
    }
    PrimitiveType const primitive_enum = PrimitiveType::_from_string_nocase(primitive_name.c_str());
    return CreatePrimitive(primitive_enum, std::forward<OptionMap>(options));
}

}  // namespace algos
