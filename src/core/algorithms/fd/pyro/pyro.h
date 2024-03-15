#pragma once

#include <list>
#include <mutex>

#include "algorithms/fd/pli_based_fd_algorithm.h"
#include "algorithms/fd/pyrocommon/core/dependency_consumer.h"
#include "algorithms/fd/pyrocommon/core/search_space.h"

namespace algos {

/* Class for mining FD with pyro algorithm */
class Pyro : public DependencyConsumer, public PliBasedFDAlgorithm {
private:
    std::list<std::unique_ptr<SearchSpace>> search_spaces_;

    CachingMethod caching_method_ = CachingMethod::kCoin;
    CacheEvictionMethod eviction_method_ = CacheEvictionMethod::kDefault;
    double caching_method_value_;

    pyro::Parameters parameters_;

    void RegisterOptions();
    void MakeExecuteOptsAvailableFDInternal() final;

    void ResetStateFd() final;
    unsigned long long ExecuteInternal() final;

public:
    Pyro();
};

}  // namespace algos
