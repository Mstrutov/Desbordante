#pragma once

class ProfilingContext;

#include "CacheEvictionMethod.h"
#include "CachingMethod.h"
#include "ProfilingContext.h"
#include "ColumnLayoutRelationData.h"

#include <mutex>

namespace util {

class PLICache {
private:
    class PositionListIndexRank {
    public:
        Vertical const* vertical_;
        std::shared_ptr<PositionListIndex> pli_;
        int added_arity_;

        PositionListIndexRank(Vertical const* vertical, std::shared_ptr<PositionListIndex> pli, int initial_arity):
            vertical_(vertical), pli_(pli), added_arity_(initial_arity) {}
    };
    //using CacheMap = VerticalMap<PositionListIndex>;
    ColumnLayoutRelationData* relation_data_;
    std::unique_ptr<VerticalMap<PositionListIndex>> index_;
    //usageCounter - for parallelism

    int saved_intersections_ = 0;

    mutable std::mutex getting_pli_mutex_;

    CachingMethod caching_method_;
    CacheEvictionMethod eviction_method_;
    double caching_method_value_;
    //long long maximumAvailableMemory_ = 0;
    double maximum_entropy_;
    double mean_entropy_;
    double min_entropy_;
    double median_entropy_;
    double median_gini_;
    double median_inverted_entropy_;

    std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> CachingProcess(Vertical const& vertical,
                                                                                        std::unique_ptr<PositionListIndex> pli,
                                                                                        ProfilingContext* profiling_context);
public:
    PLICache(ColumnLayoutRelationData* relation_data, CachingMethod caching_method, CacheEvictionMethod eviction_method,
             double caching_method_value, double min_entropy, double mean_entropy, double median_entropy,
             double maximum_entropy, double median_gini, double median_inverted_entropy);

    PositionListIndex* Get(Vertical const& vertical);
    std::variant<PositionListIndex*, std::unique_ptr<PositionListIndex>> GetOrCreateFor(
            Vertical const& vertical, ProfilingContext* profiling_context);

    void SetMaximumEntropy(double e) { maximum_entropy_ = e; }

    size_t Size() const;

    // returns ownership of single column PLIs back to ColumnLayoutRelationData
    virtual ~PLICache();
};

} // namespace util

