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
        int added_arity_;

        PositionListIndexRank(Vertical const* vertical,
                              std::shared_ptr<util::PositionListIndex> pli, int initial_arity)
            : vertical_(vertical), pli_(pli), added_arity_(initial_arity) {}
    };

    ColumnLayoutRelationData* relation_data_;
    std::unique_ptr<util::VerticalMap<util::PositionListIndex>> index_;

    int saved_intersections_ = 0;

    mutable std::mutex getting_pli_mutex_;

    CachingMethod caching_method_;
    CacheEvictionMethod eviction_method_;
    double caching_method_value_;

    double median_inverted_entropy_;

    std::variant<util::PositionListIndex*, std::unique_ptr<util::PositionListIndex>>
    CachingProcess(Vertical const& vertical, std::unique_ptr<util::PositionListIndex> pli);
public:
    PartitionStorage(ColumnLayoutRelationData* relation_data,
                     CachingMethod caching_method, CacheEvictionMethod eviction_method);

    util::PositionListIndex* Get(Vertical const& vertical);
    std::variant<util::PositionListIndex*, std::unique_ptr<util::PositionListIndex>>
    GetOrCreateFor(Vertical const& vertical);

    size_t Size() const;

    virtual ~PartitionStorage();
};
