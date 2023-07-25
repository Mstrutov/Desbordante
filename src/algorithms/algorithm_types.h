#pragma once

#include <enum.h>

#include "algorithms/algorithms.h"

namespace algos {

using AlgorithmTypes =
        std::tuple<Depminer, DFD, FastFDs, FDep, Fd_mine, Pyro, Tane, FUN, hyfd::HyFD, Aid, Apriori,
                   metric::MetricVerifier, DataStats, fd_verifier::FDVerifier, HyUCC,
                   cfd::FDFirstAlgorithm, ACAlgorithm>;

// clang-format off
/* Enumeration of all supported non-pipeline algorithms. If you implement a new
 * algorithm please add its corresponding value to this enum and to the type
 * tuple above.
 * NOTE: algorithm string name representation is taken from the value in this
 * enum, so name it appropriately (lowercase and without additional symbols).
 */
BETTER_ENUM(AlgorithmType, char,
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

/* FD verifier algorithm */
    fd_verifier,

/* Unique Column Combination mining algorithms */
    hyucc,

/* CFD mining algorithms */
    fd_first_dfs,

/* Algebraic constraints mining algorithm*/
    ac
)
// clang-format on

static_assert(std::tuple_size_v<AlgorithmTypes> == AlgorithmType::_size());

}  // namespace algos
