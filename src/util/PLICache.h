#pragma once

class ProfilingContext;

#include "CacheEvictionMethod.h"
#include "CachingMethod.h"
#include "ProfilingContext.h"
#include "ColumnLayoutRelationData.h"

#include <mutex>

class PLICache {
private:
    class PositionListIndexRank {
    public:
        Vertical const* vertical_;
        std::shared_ptr<PositionListIndex> pli_;
        int addedArity_;

        PositionListIndexRank(Vertical const* vertical, std::shared_ptr<PositionListIndex> pli, int initialArity):
            vertical_(vertical), pli_(pli), addedArity_(initialArity) {}
    };
    //using CacheMap = VerticalMap<PositionListIndex>;
    ColumnLayoutRelationData* relationData_;
    std::unique_ptr<VerticalMap<PositionListIndex>> index_;
    //usageCounter - for parallelism

    int savedIntersections_ = 0;

    mutable std::mutex gettingPLIMutex;

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

    std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> cachingProcess(Vertical const& vertical,
                        std::unique_ptr<PositionListIndex> pli,
                        ProfilingContext* profilingContext);
public:
    PLICache(ColumnLayoutRelationData* relationData, CachingMethod cachingMethod, CacheEvictionMethod evictionMethod,
             double cachingMethodValue, double minEntropy, double meanEntropy, double medianEntropy,
             double maximumEntropy, double medianGini, double medianInvertedEntropy);

    PositionListIndex* get(Vertical const& vertical);
    std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> getOrCreateFor(
            Vertical const& vertical, ProfilingContext* profilingContext);

    void setMaximumEntropy(double e) { maximumEntropy_ = e; }

    size_t size() const;

    // returns ownership of single column PLIs back to ColumnLayoutRelationData
    virtual ~PLICache();
};