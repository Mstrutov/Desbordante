#pragma once

class ProfilingContext;

#include "CacheEvictionMethod.h"
#include "CachingMethod.h"
#include "ProfilingContext.h"
#include "ColumnLayoutRelationData.h"

#include <mutex>

class PartitionStorage {
private:
    class PositionListIndexRank {
    public:
        Vertical const* vertical_;
        std::shared_ptr<PositionListIndex> pli_;
        int addedArity_;

        PositionListIndexRank(Vertical const* vertical, std::shared_ptr<PositionListIndex> pli, int initialArity):
                vertical_(vertical), pli_(pli), addedArity_(initialArity) {}
    };
    ColumnLayoutRelationData* relationData_;
    std::unique_ptr<VerticalMap<PositionListIndex>> index_;

    int savedIntersections_ = 0;

    mutable std::mutex gettingPLIMutex;

    CachingMethod cachingMethod_;
    CacheEvictionMethod evictionMethod_;
    double cachingMethodValue_;

    double medianInvertedEntropy_;

    std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> cachingProcess(Vertical const& vertical,
                                                                                        std::unique_ptr<PositionListIndex> pli);
public:
    PartitionStorage(ColumnLayoutRelationData* relationData, CachingMethod cachingMethod, CacheEvictionMethod evictionMethod);

    PositionListIndex* get(Vertical const& vertical);
    std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> getOrCreateFor(
            Vertical const& vertical);

    size_t size() const;

    virtual ~PartitionStorage();
};
