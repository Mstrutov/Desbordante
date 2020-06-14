#pragma once

#include <list>
#include <map>
#include <utility>
#include <vector>
#include "AgreeSetSample.h"

//TODO: Java long ~ C++ long long => consider replacing ints with longlongs
class ListAgreeSetSample : public AgreeSetSample {
private:
    struct Entry {

        int count_;
        std::shared_ptr<std::vector<long long>> agreeSet_;

        Entry(std::shared_ptr<std::vector<long long>> agreeSet, int count) : agreeSet_(std::move(agreeSet)), count_(count) {}
    };

    std::list<Entry> agreeSetCounters_;

    static std::shared_ptr<std::vector<long long>> bitSetToLongLongVector(boost::dynamic_bitset<> const& bitset);
public:
    static std::shared_ptr<ListAgreeSetSample> createFocusedFor(std::shared_ptr<ColumnLayoutRelationData> relation,
                                                std::shared_ptr<Vertical> restrictionVertical,
                                                std::shared_ptr<PositionListIndex> restrictionPLi,
                                                unsigned int sampleSize);

    //100% should use move semantics on agreeSetCounters.
    //in Java code relation is a reference to base class RelationData, but in fact it references to CLRD, so here latter is used
    ListAgreeSetSample(std::shared_ptr<ColumnLayoutRelationData> relation, std::shared_ptr<Vertical> focus, unsigned int sampleSize, unsigned long long populationSize,
                        std::map<boost::dynamic_bitset<>, int> const & agreeSetCounters);

    long long getNumAgreeSupersets(std::shared_ptr<Vertical> agreement) override;
    long long getNumAgreeSupersets(std::shared_ptr<Vertical> agreement, std::shared_ptr<Vertical> disagreement) override;
    std::shared_ptr<std::vector<long long>> getNumAgreeSupersetsExt(std::shared_ptr<Vertical> agreement, std::shared_ptr<Vertical> disagreement) override;
};