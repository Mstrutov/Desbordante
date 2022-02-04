#pragma once

#include <map>
#include <utility>
#include <vector>
#include "AgreeSetSample.h"

namespace util {

//TODO: Java long ~ C++ long long => consider replacing ints with longlongs
class ListAgreeSetSample : public AgreeSetSample {
private:
    struct Entry {

        unsigned int count;
        std::shared_ptr<std::vector<unsigned long long>> agree_set;

        Entry(std::shared_ptr<std::vector<unsigned long long>> agree_set, unsigned int count)
        : count(count), agree_set(std::move(agree_set)) {}
    };

    std::vector<Entry> agree_set_counters_;

public:
    static std::unique_ptr<std::vector<unsigned long long>> BitSetToLongLongVector(
            boost::dynamic_bitset<> const& bitset);

    static std::unique_ptr<ListAgreeSetSample> CreateFocusedFor(
        ColumnLayoutRelationData const* relation, Vertical const& restriction_vertical,
        PositionListIndex const* restriction_p_li, unsigned int sample_size, CustomRandom& random);

    ListAgreeSetSample(ColumnLayoutRelationData const* relation, Vertical const& focus,
                       unsigned int sample_size, unsigned long long population_size,
                       std::unordered_map<boost::dynamic_bitset<>, int> const& agree_set_counters);

    unsigned long long GetNumAgreeSupersets(Vertical const& agreement) const override;
    unsigned long long GetNumAgreeSupersets(Vertical const& agreement, Vertical const& disagreement) const override;
    std::unique_ptr<std::vector<unsigned long long>> GetNumAgreeSupersetsExt(
            Vertical const& agreement, Vertical const& disagreement) const override;
};

} // namespace util

