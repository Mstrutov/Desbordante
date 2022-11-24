#pragma once

#include <boost/mp11.hpp>
#include <enum.h>

#include "algorithms/algorithms.h"

namespace algos {

using PrimitiveTypes = std::tuple<Depminer, DFD, FastFDs, FDep, Fd_mine, Pyro, Tane, FUN,
        hyfd::HyFD, Aid, Apriori, metric::MetricVerifier, DataStats, CFDDiscovery>;

/* Enumeration of all supported non-pipeline primitives. If you implement a new
 * primitive please add its corresponding value to this enum and to the type
 * tuple above.
 * NOTE: algorithm string name representation is taken from the value in this
 * enum, so name it appropriately (lowercase and without additional symbols).
 */
BETTER_ENUM(PrimitiveType, char,
/* Functional dependency mining algorithms */
    depminer = 0,
    dfd,
    fastfds,
    fdep,
    fdmine,
    pyro,
    tane,
    fun,
    hyfd,
    aidfd,

/* Association rules mining algorithms */
    apriori,

/* Metric verifier algorithm */
    metric,

/* Statistic algorithms */
    stats,

/* Conditional functional dependencies mining algorithm */
    cfd
)

template <typename PrimitiveBase = Primitive>
std::unique_ptr<PrimitiveBase> CreatePrimitiveInstance(PrimitiveType primitive) {
    auto const create = [](auto I) -> std::unique_ptr<PrimitiveBase> {
        using PrimitiveType = std::tuple_element_t<I, PrimitiveTypes>;
        if constexpr (std::is_convertible_v<PrimitiveType *, PrimitiveBase *>) {
            return std::make_unique<PrimitiveType>();
        }
        else {
            throw std::invalid_argument("Cannot use "
                                        + boost::typeindex::type_id<PrimitiveType>().pretty_name()
                                        + " as "
                                        + boost::typeindex::type_id<PrimitiveBase>().pretty_name());
        }
    };

    return boost::mp11::mp_with_index<std::tuple_size<PrimitiveTypes>>(
            static_cast<size_t>(primitive), create);
}

template <typename PrimitiveBase>
std::vector<PrimitiveType> GetAllDerived() {
    auto const is_derived = [](auto I) -> bool {
        using PrimitiveType = std::tuple_element_t<I, PrimitiveTypes>;
        return std::is_base_of_v<PrimitiveBase, PrimitiveType>;
    };
    std::vector<PrimitiveType> derived_from_base{};
    for (PrimitiveType prim : PrimitiveType::_values()) {
        if (boost::mp11::mp_with_index<std::tuple_size<PrimitiveTypes>>(static_cast<size_t>(prim),
                                                                        is_derived)) {
            derived_from_base.push_back(prim);
        }
    }
    return derived_from_base;
}

}  // namespace algos
