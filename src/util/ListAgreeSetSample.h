#pragma once

#include <map>
#include <utility>
#include <vector>
#include "AgreeSetSample.h"

//TODO: Java long ~ C++ long long => consider replacing ints with longlongs
class ListAgreeSetSample : public AgreeSetSample {
private:
    struct Entry {

        unsigned int count_;
        std::shared_ptr<std::vector<unsigned long long>> agreeSet_;

        Entry(std::shared_ptr<std::vector<unsigned long long>> agreeSet, unsigned int count)
        : count_(count), agreeSet_(std::move(agreeSet)) {}
    };

    std::vector<Entry> agreeSetCounters_;

public:
    static std::unique_ptr<std::vector<unsigned long long>> bitSetToLongLongVector(
            boost::dynamic_bitset<> const& bitset);

    static std::unique_ptr<ListAgreeSetSample> createFocusedFor(
            ColumnLayoutRelationData* relation, Vertical const& restrictionVertical, PositionListIndex* restrictionPLi,
            unsigned int sampleSize, CustomRandom& random);

    ListAgreeSetSample(ColumnLayoutRelationData* relation, Vertical const& focus,
                       unsigned int sampleSize, unsigned long long populationSize,
                       std::unordered_map<boost::dynamic_bitset<>, int> const& agreeSetCounters);

    unsigned long long getNumAgreeSupersets(Vertical const& agreement) const override;
    unsigned long long getNumAgreeSupersets(Vertical const& agreement, Vertical const& disagreement) const override;
    std::unique_ptr<std::vector<unsigned long long>> getNumAgreeSupersetsExt(
            Vertical const& agreement, Vertical const& disagreement) const override;
};