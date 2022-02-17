#pragma once

#include <list>
#include <mutex>

#include "CSVParser.h"
#include "DependencyConsumer.h"
#include "PliBasedFDAlgorithm.h"
#include "SearchSpace.h"

class Pyro : public DependencyConsumer, public PliBasedFDAlgorithm {
private:
    constexpr static const char* kSeed = "seed";
    constexpr static const char* kMaxError = "error";

    std::list<std::unique_ptr<SearchSpace>> search_spaces_;

    CachingMethod caching_method_ = CachingMethod::kCoin;
    CacheEvictionMethod eviction_method_ = CacheEvictionMethod::kDefault;
    double caching_method_value_;

    Configuration configuration_;

    unsigned long long ExecuteInternal() override;
    void init();

public:
    explicit Pyro(Config const& config);
    explicit Pyro(std::shared_ptr<ColumnLayoutRelationData> relation, Config const& config);
};
