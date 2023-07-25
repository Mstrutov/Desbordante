#pragma once

#include <list>
#include <mutex>

#include "algorithms/pli_based_fd_algorithm.h"
#include "core/dependency_consumer.h"
#include "core/search_space.h"

namespace algos {

class Pyro : public DependencyConsumer, public PliBasedFDAlgorithm {
private:
    std::list<std::unique_ptr<SearchSpace>> search_spaces_;

    CachingMethod caching_method_ = CachingMethod::kCoin;
    CacheEvictionMethod eviction_method_ = CacheEvictionMethod::kDefault;
    double caching_method_value_;

    Parameters parameters_;

    void RegisterOptions();

    void ResetStateFd() final;
    unsigned long long ExecuteInternal() final;

public:
    Pyro(bool request_prepared_data = false);
};

}  // namespace algos
