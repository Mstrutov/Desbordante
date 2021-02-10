#pragma once

class ProfilingContext;
#include "CacheEvictionMethod.h"
#include "CachingMethod.h"
#include "ProfilingContext.h"
#include "ColumnLayoutRelationData.h"

class PLICache {
private:
    class PositionListIndexRank {
    public:
        std::shared_ptr<Vertical> vertical_;
        std::shared_ptr<PositionListIndex> pli_;
        int addedArity_;

        PositionListIndexRank(std::shared_ptr<Vertical> vertical, std::shared_ptr<PositionListIndex> pli, int initialArity):
            vertical_(vertical), pli_(pli), addedArity_(initialArity) {}
    };
    using CacheMap = VerticalMap<std::shared_ptr<PositionListIndex>>;
    std::weak_ptr<ColumnLayoutRelationData> relationData_;
    std::shared_ptr<CacheMap> index_;    //unique_ptr?
    //usageCounter - for parallelism

    int savedIntersections_ = 0;

    CachingMethod cachingMethod_;
    CacheEvictionMethod evictionMethod_;
    double cachingMethodValue_;
    //long long maximumAvailableMemory_ = 0;
    double maximumEntropy_;
    double meanEntropy_;
    double minEntropy_;
    double medianEntropy_;
    double medianGini_;
    double medianInvertedEntropy_;

    void cachingProcess(Vertical const& vertical, std::shared_ptr<PositionListIndex> pli, ProfilingContext* profilingContext);
public:
    PLICache(std::shared_ptr<ColumnLayoutRelationData> relationData, CachingMethod cachingMethod, CacheEvictionMethod evictionMethod,
             double cachingMethodValue, double minEntropy, double meanEntropy, double medianEntropy, double maximumEntropy, double medianGini, double medianInvertedEntropy);

    std::shared_ptr<PositionListIndex> get(Vertical const& vertical);
    std::shared_ptr<PositionListIndex> getOrCreateFor(Vertical const& vertical, ProfilingContext* profilingContext);

    void setMaximumEntropy(double e) { maximumEntropy_ = e; }

    size_t size() const;
};