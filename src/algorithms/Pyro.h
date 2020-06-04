#pragma once
#include <list>
#include "DependencyConsumer.h"
#include "SearchSpace.h"
#include "CSVParser.h"

class Pyro : public DependencyConsumer {
private:
    CSVParser inputGenerator_;
    std::list<std::shared_ptr<SearchSpace>> searchSpaces_;
    //config, uccConsumer, fdConsumer,

    CachingMethod cachingMethod_;
    CacheEvictionMethod evictionMethod_;
    double cachingMethodValue;

    Configuration configuration_;
public:
    explicit Pyro(fs::path const& path);
    void execute();
};