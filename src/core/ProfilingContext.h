#pragma once

#include <random>
//forward declaration
class PLICache;

template<class Value>
class VerticalMap;
#include "core/Configuration.h"
#include "util/AgreeSetSample.h"
//#include "util/PLICache.h"


//Dependency Consumer?
class ProfilingContext {
public:
    enum class ObjectToCache {
        PLI,
        AS
    };

    Configuration configuration_;
    std::shared_ptr<PLICache> pliCache_;            //unique_ptr?
    std::shared_ptr<VerticalMap<std::shared_ptr<AgreeSetSample>>> agreeSetSamples_;     //unique_ptr?

    //std::mt19937 random_;
};