#pragma once

#include <boost/mp11.hpp>

#include "algorithms/algorithm_types.h"

namespace algos {

template <typename AlgorithmBase = Algorithm, typename... ConstructorArgs>
std::unique_ptr<AlgorithmBase> CreateAlgorithmInstance(AlgorithmType algorithm,
                                                       ConstructorArgs&&... args) {
    auto const create = [&args...](auto I) -> std::unique_ptr<AlgorithmBase> {
        using AlgorithmType = std::tuple_element_t<I, AlgorithmTypes>;
        if constexpr (std::is_convertible_v<AlgorithmType*, AlgorithmBase*>) {
            return std::make_unique<AlgorithmType>(std::forward<ConstructorArgs>(args)...);
        } else {
            throw std::invalid_argument(
                    "Cannot use " + boost::typeindex::type_id<AlgorithmType>().pretty_name() +
                    " as " + boost::typeindex::type_id<AlgorithmBase>().pretty_name());
        }
    };

    return boost::mp11::mp_with_index<std::tuple_size<AlgorithmTypes>>(
            static_cast<size_t>(algorithm), create);
}

template <typename Base>
bool IsBaseOf(AlgorithmType algorithm) {
    auto const is_derived = [](auto I) -> bool {
        using AlgoType = std::tuple_element_t<I, AlgorithmTypes>;
        return std::is_base_of_v<Base, AlgoType>;
    };
    return boost::mp11::mp_with_index<std::tuple_size<AlgorithmTypes>>(
            static_cast<size_t>(algorithm), is_derived);
}

template <typename AlgorithmBase>
std::vector<AlgorithmType> GetAllDerived() {
    std::vector<AlgorithmType> derived_from_base{};
    for (AlgorithmType algo : AlgorithmType::_values()) {
        if (IsBaseOf<AlgorithmBase>(algo)) {
            derived_from_base.push_back(algo);
        }
    }
    return derived_from_base;
}

}  // namespace algos
