#pragma once
#include <list>
#include <mutex>

#include "CSVParser.h"
#include "DependencyConsumer.h"
#include "PliBasedFDAlgorithm.h"
#include "SearchSpace.h"

class Pyro : public DependencyConsumer, public PliBasedFDAlgorithm {
private:
    std::list<std::unique_ptr<SearchSpace>> searchSpaces_;

    CachingMethod cachingMethod_;
    CacheEvictionMethod evictionMethod_;
    double cachingMethodValue;

    Configuration configuration_;

    unsigned long long executeInternal() override;
public:
    explicit Pyro(std::filesystem::path const& path, char separator = ',', bool hasHeader = true,
                  int seed = 0, double maxError = 0, unsigned int maxLHS = -1, int parallelism = 0);
};
