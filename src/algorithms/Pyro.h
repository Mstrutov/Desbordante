#pragma once
#include <list>
#include <mutex>

#include "CSVParser.h"
#include "DependencyConsumer.h"
#include "PliBasedFDAlgorithm.h"
#include "SearchSpace.h"

class Pyro : public DependencyConsumer, public PliBasedFDAlgorithm {
private:
    std::list<std::unique_ptr<SearchSpace>> search_spaces_;

    CachingMethod caching_method_;
    CacheEvictionMethod eviction_method_;
    double caching_method_value_;

    Configuration configuration_;

    unsigned long long ExecuteInternal() override;
public:
    explicit Pyro(std::filesystem::path const& path, char separator = ',', bool has_header = true,
                  int seed = 0, double max_error = 0, unsigned int max_lhs = -1,
                  int parallelism = 0);
};
