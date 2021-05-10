#pragma once

#include <set>
#include <deque>
#include <vector>

#include "Vertical.h"
#include "ColumnLayoutRelationData.h"

using AgreeSet = Vertical;

enum class AgreeSetsGenMethod {
    kUsingVectorOfIDSets = 0,
    kUsingMapOfIDSets       ,
    kUsingGetAgreeSet       , // Simplest method by Wyss et al
    kUsingMCAndGetAgreeSet  , // From maximal representation using getAgreeSet()
};

class AgreeSetFactory {
public:
    explicit AgreeSetFactory(ColumnLayoutRelationData const* const rel)
        : relation_(rel) {}

    ColumnLayoutRelationData const* getRelation() const { return relation_; }

    // Computes all agree sets of `relation_` using specified method
    template<AgreeSetsGenMethod method = AgreeSetsGenMethod::kUsingVectorOfIDSets>
    std::set<AgreeSet> genAgreeSets() const;
    std::set<std::vector<int>> genPLIMaxRepresentation() const;
    AgreeSet getAgreeSet(int const tuple1_index, int const tuple2_index) const;

private:
    void calculateSupersets(std::set<std::vector<int>>& max_representation,
                            std::deque<std::vector<int>> partition) const;
    ColumnLayoutRelationData const* const relation_;
};