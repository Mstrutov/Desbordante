#pragma once

class ProfilingContext;

#include "CacheEvictionMethod.h"
#include "CachingMethod.h"
#include "ProfilingContext.h"
#include "ColumnLayoutRelationData.h"
#include "util/VerticalMap.h"

#include <mutex>

class PartitionStorage {
private:
    class PositionListIndexRank {
    public:
        Vertical const* vertical_;
        std::shared_ptr<util::PositionListIndex> pli_;
        int addedArity_;

        PositionListIndexRank(Vertical const* vertical,
                              std::shared_ptr<util::PositionListIndex> pli,
                              int initialArity):
                vertical_(vertical), pli_(pli), addedArity_(initialArity) {}
    };
    ColumnLayoutRelationData* relationData_;
    std::unique_ptr<util::VerticalMap<util::PositionListIndex>> index_;

    int savedIntersections_ = 0;

    mutable std::mutex gettingPLIMutex;

    CachingMethod cachingMethod_;
    CacheEvictionMethod evictionMethod_;
    double cachingMethodValue_;

    double medianInvertedEntropy_;

    std::variant<util::PositionListIndex*, std::unique_ptr<util::PositionListIndex>>
    cachingProcess(Vertical const& vertical, std::unique_ptr<util::PositionListIndex> pli);
public:
    PartitionStorage(ColumnLayoutRelationData* relationData,
                     CachingMethod cachingMethod, CacheEvictionMethod evictionMethod);

    util::PositionListIndex* get(Vertical const& vertical);
    std::variant<util::PositionListIndex*, std::unique_ptr<util::PositionListIndex>>
    getOrCreateFor(Vertical const& vertical);

    size_t size() const;

    virtual ~PartitionStorage();
};
