#pragma once

#include <iostream>
#include <random>
#include <string>

#include "AgreeSetSample.h"
#include "CacheEvictionMethod.h"
#include "CachingMethod.h"
#include "Configuration.h"
#include "custom/CustomRandom.h"
#include "PartialFD.h"
#include "PartialKey.h"
#include "DependencyConsumer.h"

//forward declaration
class PLICache;

template<class Value>
class VerticalMap;

//Dependency Consumer?
class ProfilingContext : public DependencyConsumer {
private:
    Configuration configuration_;
    std::unique_ptr<PLICache> pliCache_;
    std::unique_ptr<VerticalMap<AgreeSetSample>> agreeSetSamples_;     //unique_ptr?
    ColumnLayoutRelationData* relationData_;
    std::mt19937 random_;
    CustomRandom customRandom_;

    AgreeSetSample const* createColumnFocusedSample(
            Vertical const& focus, PositionListIndex const* restrictionPLI, double boostFactor);

public:
    enum class ObjectToCache {
        PLI,
        AS
    };

    ProfilingContext(Configuration  configuration, ColumnLayoutRelationData* relationData,
            std::function<void (PartialKey const&)> const& uccConsumer,
            std::function<void (PartialFD const&)> const& fdConsumer,
            CachingMethod const& cachingMethod, CacheEvictionMethod const& evictionMethod, double cachingMethodValue);

    // Non-const as RandomGenerator state gets changed
    AgreeSetSample const* createFocusedSample(Vertical const& focus, double boostFactor);
    std::shared_ptr<AgreeSetSample const> getAgreeSetSample(Vertical const& focus) const;
    PLICache* getPLICache() { return pliCache_.get(); }
    bool isAgreeSetSamplesEmpty() const { return agreeSetSamples_ == nullptr; }
    RelationalSchema const* getSchema() const { return relationData_->getSchema(); }

    Configuration const& getConfiguration() const { return configuration_; }
    ColumnLayoutRelationData const* getColumnLayoutRelationData() const { return relationData_; }
    PLICache const* getPLICache() const { return pliCache_.get(); }

    // get int in range [0, upperBound) from the uniform distribution
    // int nextInt(int upperBound) { return std::uniform_int_distribution<int>{0, upperBound}(random_); }
    int nextInt(int upperBound) { return customRandom_.nextInt(upperBound); }
    double nextDouble() { return customRandom_.nextDouble(); }

    ~ProfilingContext();

    static double getMaximumEntropy(ColumnLayoutRelationData const* cd1);
    static double getMinEntropy(ColumnLayoutRelationData const* cd1);
    static double getMedianEntropy(ColumnLayoutRelationData const* relationData);
    static double getMedianInvertedEntropy(ColumnLayoutRelationData const* relationData);
    static double getMeanEntropy(ColumnLayoutRelationData const* relationData);
    static double getMedianGini(ColumnLayoutRelationData const* relationData);
private:
    static double getMedianValue(std::vector<double> && values, std::string const& measureName);
    static double setMaximumEntropy(ColumnLayoutRelationData const* relationData, CachingMethod const & cachingMethod);
};
