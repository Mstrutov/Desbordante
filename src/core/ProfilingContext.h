#pragma once

#include <random>
//forward declaration
class PLICache;

template<class Value>
class VerticalMap;
#include "core/Configuration.h"
#include "util/AgreeSetSample.h"
#include "model/PartialFD.h"
#include "model/PartialKey.h"
#include "caching/CacheEvictionMethod.h"
#include "caching/CachingMethod.h"
#include "util/AgreeSetSample.h"
#include "core/DependencyConsumer.h"
//#include "util/PLICache.h"


//Dependency Consumer?
class ProfilingContext : public DependencyConsumer {
public:
    enum class ObjectToCache {
        PLI,
        AS
    };

    Configuration configuration_;
    std::shared_ptr<PLICache> pliCache_;            //unique_ptr?
    std::shared_ptr<VerticalMap<std::shared_ptr<AgreeSetSample>>> agreeSetSamples_;     //unique_ptr?

    std::shared_ptr<ColumnLayoutRelationData> relationData_;
    std::mt19937 random_;

    //std::function<void (PartialFD const&)> fdConsumer_;
    //std::function<void (PartialKey const&)> uccConsumer_;

    // initialize random_ using std::random_device
    ProfilingContext(Configuration const& configuration, std::shared_ptr<ColumnLayoutRelationData> relationData,
            std::function<void (PartialKey const&)> const& uccConsumer, std::function<void (PartialFD const&)> const& fdConsumer,
            CachingMethod const& cachingMethod, CacheEvictionMethod const& evictionMethod, double cachingMethodValue);

    std::shared_ptr<AgreeSetSample> createFocusedSample(std::shared_ptr<Vertical> focus, double boostFactor);
    // Retrieve an AgreeSetSample with a best possible sampling ratio
    std::shared_ptr<AgreeSetSample> getAgreeSetSample(std::shared_ptr<Vertical> focus);
    std::shared_ptr<RelationalSchema> getSchema() { return relationData_->getSchema(); }

    // get int in range [0, upperBound) from the uniform distribution
    int nextInt(int upperBound) { return std::uniform_int_distribution<int>{0, upperBound}(random_); }

    static double getMaximumEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData);
    static double getMinEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData);
    static double getMedianEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData);
    static double getMedianInvertedEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData);
    static double getMeanEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData);
    static double getMedianGini(std::shared_ptr<ColumnLayoutRelationData> relationData);
private:
    static double setMaximumEntropy(std::shared_ptr<ColumnLayoutRelationData> relationData, CachingMethod const & cachingMethod);
};