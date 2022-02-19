#pragma once

#include <boost/any.hpp>

#include "better_enums/enum.h"

#include "Algorithms.h"
#include "TypoMiner.h"

namespace algos {

BETTER_ENUM(AlgoMiningType, char,
#if 1
    fd = 0,
    typos
#else
    fd = 0, /* Functional dependency mining */
    cfd,    /* Conditional functional dependency mining */
    ar,     /* Association rule mining */
    key,    /* Key mining */
    typos   /* Typo mining */
#endif
);

/* Enumeration of all supported algorithms. If you implemented new algorithm
 * please add new corresponding value to this enum.
 * NOTE: algorithm string name represenation is taken from value in this enum,
 * so name it appropriately (lowercase and without additional symbols).
 */
BETTER_ENUM(Algo, char,
    /* Functional dependency mining algorithms */
    depminer = 0,
    dfd,
    fastfds,
    fdep,
    fdmine,
    pyro,
    tane
);

using StdParamsMap = std::unordered_map<std::string, boost::any>;
using AlgorithmTypesTuple = std::tuple<Depminer, DFD, FastFDs, FDep, Fd_mine, Pyro, Tane>;

namespace details {

/* Consider using boost:hana here to get rid of this nightmare */
template <size_t... Is, typename F>
constexpr void StaticForImpl(F&& f, std::index_sequence<Is...>) {
    (f(std::integral_constant<size_t, Is>()), ...);
}

template <size_t N, typename F>
constexpr void StaticFor(F&& f) {
    StaticForImpl(f, std::make_index_sequence<N>());
}

template <typename Tuple, typename F>
constexpr void ForAllTypes(F&& f) {
    StaticFor<std::tuple_size_v<Tuple>>([&f](auto N) {
        using TupleElement = std::tuple_element_t<N, Tuple>;

        f((TupleElement const*)nullptr, N);
    });
}

template <typename R, typename Tuple, typename F>
constexpr R SelectType(F&& f, size_t i) {
    R r;
    bool found = false;

    ForAllTypes<Tuple>([&r, &found, &f, i](auto const* p, auto N) {
        if (i == N) {
            r = f(p);
            found = true;
        }
    });

    if (!found) {
        throw std::invalid_argument("Wrong enum value");
    }

    return r;
}

template <typename AlgorithmBase = Primitive, typename AlgorithmsToSelect = AlgorithmTypesTuple,
          typename EnumType, typename... Args>
auto CreatePrimitiveInstanceImpl(EnumType const enum_value, Args&&... args) {
    return SelectType<std::unique_ptr<AlgorithmBase>, AlgorithmsToSelect>(
        [&args...](auto const* p) {
            using AlgorithmType = std::decay_t<decltype(*p)>;
            return std::make_unique<AlgorithmType>(std::forward<Args>(args)...);
        },
        static_cast<size_t>(enum_value));
}

template <template <typename...> class Wrapper, typename AlgorithmBase = Primitive,
          typename AlgorithmsToWrap = AlgorithmTypesTuple, typename EnumType, typename... Args>
auto CreateAlgoWrapperInstanceImpl(EnumType const enum_value, Args&&... args) {
    static_assert(std::is_base_of_v<AlgorithmBase, Wrapper<AlgorithmBase>>,
                  "Wrapper should be derived from AlgorithmBase");
    return SelectType<std::unique_ptr<AlgorithmBase>, AlgorithmsToWrap>(
        [&args...](auto const* p) {
            using AlgorithmType = std::decay_t<decltype(*p)>;
            return std::make_unique<Wrapper<AlgorithmType>>(std::forward<Args>(args)...);
        },
        static_cast<size_t>(enum_value));
}

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
std::unique_ptr<Primitive> CreateFDAlgorithmInstance(Algo const algo, ParamsMap&& params) {
    FDAlgorithm::Config const config =
        CreateFDAlgorithmConfigFromMap(std::forward<ParamsMap>(params));

    return details::CreatePrimitiveInstanceImpl(algo, config);
}

template <typename ParamsMap>
std::unique_ptr<Primitive> CreateTyposMinerInstance(Algo const algo, ParamsMap&& params) {
    /* Typos miner has FDAlgorithm configuration */
    FDAlgorithm::Config const config =
        CreateFDAlgorithmConfigFromMap(std::forward<ParamsMap>(params));

    return details::CreateAlgoWrapperInstanceImpl<TypoMiner>(algo, config);
}

} // namespace details

template <typename ParamsMap>
std::unique_ptr<Primitive> CreateAlgorithmInstance(AlgoMiningType const task, Algo const algo,
                                                   ParamsMap&& params) {
    switch (task) {
    case AlgoMiningType::fd:
        return details::CreateFDAlgorithmInstance(algo, std::forward<ParamsMap>(params));
    case AlgoMiningType::typos:
        return details::CreateTyposMinerInstance(algo, std::forward<ParamsMap>(params));
    default:
        throw std::logic_error(task._to_string() + std::string(" task type is not supported yet."));
    }
}

template <typename ParamsMap>
std::unique_ptr<Primitive> CreateAlgorithmInstance(std::string const& task_name,
                                                   std::string const& algo_name,
                                                   ParamsMap&& params) {
    AlgoMiningType const task = AlgoMiningType::_from_string(task_name.c_str());
    Algo const algo = Algo::_from_string(algo_name.c_str());
    return CreateAlgorithmInstance(task, algo, std::forward<ParamsMap>(params));
}

}  // namespace algos
