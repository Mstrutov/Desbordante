#pragma once
#include <list>
#include "core/DependencyConsumer.h"
#include "core/SearchSpace.h"
#include "parser/CSVParser.h"

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