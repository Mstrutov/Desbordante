#pragma once

#include "better_enums/enum.h"

#include "algorithms.h"

namespace algos {

BETTER_ENUM(AlgoMiningType, char,
    kNotSpecified = -1,

    fd,     /* Functional dependency mining */
    cfd,    /* Conditional functional dependency mining */
    ar,     /* Association rule mining */
    key,    /* Keys mining */
    error   /* Errors mining */
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

std::unique_ptr<FDAlgorithm> CreateFDAlgorithmInstance(Algo const algo,
                                                       std::string const& dataset,
                                                       char const sep,
                                                       bool const has_header,
                                                       double const error,
                                                       unsigned int const max_lhs,
                                                       unsigned int const parallelism,
                                                       int seed);

/* TODO(polyntsov): generic function to create algorithm of every AlgoMiningType type.
 * Call apropriate `Create'Type'AlgorithmInstance` function?
 */
template<typename... Params>
std::unique_ptr<FDAlgorithm> CreateAlgorithmInstance(std::string const& name,
                                                     Params&&... params) {
    Algo const algo = Algo::_from_string(name.c_str());
    return CreateFDAlgorithmInstance(algo, std::forward<Params>(params)...);
}

} // namespace algos

