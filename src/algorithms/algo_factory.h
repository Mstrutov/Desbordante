#pragma once

#include <boost/any.hpp>

#include "better_enums/enum.h"

#include "algorithms.h"

namespace algos {

BETTER_ENUM(AlgoMiningType, char,
#if 1
    fd
#else
    fd,     /* Functional dependency mining */
    cfd,    /* Conditional functional dependency mining */
    ar,     /* Association rule mining */
    key,    /* Keys mining */
    error   /* Errors mining */
#endif
);

/* Enumeration of all supported algorithms. If you implemented new algorithm
 * please add new corresponding value to this enum.
 * NOTE: algorithm string name represenation is taken from value in this enum,
 * so name it appropriately (lowercase and without additional symbols).
 */
BETTER_ENUM(Algo, char,
    /* Functional dependency mining algorithms */
    depminer,
    dfd,
    fastfds,
    fdep,
    fdmine,
    pyro,
    tane
);

using StdParamsMap = std::unordered_map<std::string, boost::any>;

template <typename T, typename ParamsMap>
T ExtractParamFromMap(ParamsMap& params, std::string const& param_name) {
    auto it = params.find(param_name);
    if (it == params.end()) {
        /* Throw an exception here? Or validate parameters somewhere else (main?)? */
        assert(0);
    }
    if constexpr (std::is_same_v<ParamsMap, StdParamsMap>) {
        return boost::any_cast<T>(params.extract(it).mapped());
    } else {
        return params.extract(it).mapped().template as<T>();
    }
}

/* Really cumbersome, also copying parameter names and types throughout the project
 * (here, main, pyro and tane params). It will be really hard to maintain, need to fix this
 */
template <typename ParamsMap>
FDAlgorithm::Config CreateFDAlgorithmConfigFromMap(ParamsMap params) {
    FDAlgorithm::Config c;

    c.data = std::filesystem::current_path() / "inputData" /
             ExtractParamFromMap<std::string>(params, "data");
    c.separator = ExtractParamFromMap<char>(params, "separator");
    c.has_header = ExtractParamFromMap<bool>(params, "has_header");
    c.is_null_equal_null = ExtractParamFromMap<bool>(params, "is_null_equal_null");
    c.max_lhs = ExtractParamFromMap<unsigned int>(params, "max_lhs");
    c.parallelism = ExtractParamFromMap<ushort>(params, "threads");

    /* Is it correct to insert all specified parameters into the algorithm config, and not just the
     * necessary ones? It is definitely simpler, so for now leaving it like this
     */
    for (auto it = params.begin(); it != params.end();) {
        auto node = params.extract(it++);
        if constexpr (std::is_same_v<ParamsMap, StdParamsMap>) {
            c.special_params.emplace(std::move(node.key()), std::move(node.mapped()));
        } else {
            /* ParamsMap == boost::program_options::variable_map */
            c.special_params.emplace(std::move(node.key()), std::move(node.mapped().value()));
        }
    }

    return c;
}

template <typename ParamsMap>
std::unique_ptr<FDAlgorithm> CreateFDAlgorithmInstance(Algo const algo, ParamsMap const& params) {
    FDAlgorithm::Config const config = CreateFDAlgorithmConfigFromMap(params);

    /* A little bit ugly, but I don't know how to do it better */
    switch (algo) {
    case Algo::depminer:
        return std::make_unique<Depminer>(config);
    case Algo::dfd:
        return std::make_unique<DFD>(config);
    case Algo::fastfds:
        return std::make_unique<FastFDs>(config);
    case Algo::fdep:
        return std::make_unique<FDep>(config);
    case Algo::fdmine:
        return std::make_unique<Fd_mine>(config);
    case Algo::pyro:
        return std::make_unique<Pyro>(config);
    case Algo::tane:
        return std::make_unique<Tane>(config);
    default:
        /* Unreachable code */
        assert(0);
        break;
    }

    return nullptr;
}

template <typename ParamsMap>
std::unique_ptr<FDAlgorithm> CreateAlgorithmInstance(AlgoMiningType const task, Algo const algo,
                                                     ParamsMap const& params) {
    switch (task) {
    case AlgoMiningType::fd:
        return CreateFDAlgorithmInstance(algo, params);
        break;
    default:
        throw std::runtime_error(task._to_string() +
                                 std::string(" task type is not supported yet."));
    }
}

/* TODO(polyntsov): generic function to create algorithm of every AlgoMiningType type.
 * Call apropriate `Create'Type'AlgorithmInstance` function?
 */
template <typename ParamsMap>
std::unique_ptr<FDAlgorithm> CreateAlgorithmInstance(std::string const& task_name,
                                                     std::string const& algo_name,
                                                     ParamsMap const& params) {
    AlgoMiningType const task = AlgoMiningType::_from_string(task_name.c_str());
    Algo const algo = Algo::_from_string(algo_name.c_str());
    return CreateAlgorithmInstance(task, algo, params);
}

} // namespace algos

