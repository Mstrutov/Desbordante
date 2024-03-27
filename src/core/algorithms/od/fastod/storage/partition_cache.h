#pragma once

#include <memory>

#include "algorithms/od/fastod/model/attribute_set.h"
#include "algorithms/od/fastod/partitions/complex_stripped_partition.h"
#include "cache_with_limit.h"

namespace algos::fastod {

class PartitionCache {
private:
    CacheWithLimit<AttributeSet, ComplexStrippedPartition> cache_{static_cast<size_t>(1e8)};

    void CallProductWithAttribute(ComplexStrippedPartition& partition, size_t attribute) {
        partition.Product(attribute);

        if (partition.ShouldBeConvertedToStrippedPartition()) {
            partition.ToStrippedPartition();
        }
    }

    bool CallProductWithAttributesInCache(ComplexStrippedPartition& result,
                                          AttributeSet const& attribute_set) {
        bool is_product_called = false;

        attribute_set.Iterate(
                [this, &attribute_set, &result, &is_product_called](model::ColumnIndex attr) {
                    AttributeSet one_less = DeleteAttribute(attribute_set, attr);

                    if (one_less.Any() && cache_.Contains(one_less)) {
                        result = cache_.Get(one_less);
                        CallProductWithAttribute(result, attr);
                        is_product_called = true;
                    }
                });

        return is_product_called;
    }

public:
    void Clear() {
        cache_.Clear();
    }

    ComplexStrippedPartition GetStrippedPartition(AttributeSet const& attribute_set,
                                                  std::shared_ptr<DataFrame> data) {
        if (cache_.Contains(attribute_set)) {
            return cache_.Get(attribute_set);
        }

        ComplexStrippedPartition result_partition;
        bool is_product_called = CallProductWithAttributesInCache(result_partition, attribute_set);

        if (!is_product_called) {
            result_partition = data->IsAttributesMostlyRangeBased(attribute_set)
                                       ? ComplexStrippedPartition::Create<true>(data)
                                       : ComplexStrippedPartition::Create<false>(data);

            attribute_set.Iterate([this, &result_partition](model::ColumnIndex attr) {
                CallProductWithAttribute(result_partition, attr);
            });
        }

        cache_.Set(attribute_set, result_partition);
        return result_partition;
    }
};

}  // namespace algos::fastod
