//
// Created by alex on 22.04.2021.
//

#pragma once

#include "CacheEvictionMethod.h"
#include "CachingMethod.h"
#include "ProfilingContext.h"
#include "ColumnLayoutRelationData.h"

class PartitionStorage {
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
    //long long maximumAvailableMemory_ = 0;

    //void cachingProcess(Vertical const& vertical, std::shared_ptr<PositionListIndex> pli);

public:
    PartitionStorage(std::shared_ptr<ColumnLayoutRelationData> relationData, CachingMethod cachingMethod, CacheEvictionMethod evictionMethod);

    std::shared_ptr<PositionListIndex> get(Vertical const& vertical);
    std::shared_ptr<PositionListIndex> getOrCreateFor(Vertical const& vertical);

    size_t size() const;
};
