#pragma once

#include <set>
#include <deque>
#include <vector>
#include <unordered_set>

#include <boost/functional/hash.hpp>

#include "Vertical.h"
#include "ColumnLayoutRelationData.h"
#include "custom/CustomHashes.h"

using AgreeSet = Vertical;

enum class AgreeSetsGenMethod {
    kUsingVectorOfIDSets = 0,
    kUsingMapOfIDSets       ,
    kUsingGetAgreeSet       , // Simplest method by Wyss et al
    kUsingMCAndGetAgreeSet  , // From maximal representation using getAgreeSet()
};

class AgreeSetFactory {
public:
    using SetOfVectors = std::unordered_set<std::vector<int>, boost::hash<std::vector<int>>>;
    using SetOfAgreeSets = std::unordered_set<AgreeSet>;

    explicit AgreeSetFactory(ColumnLayoutRelationData const* const rel)
        : relation_(rel) {}

    ColumnLayoutRelationData const* getRelation() const { return relation_; }

    // Computes all agree sets of `relation_` using specified method
    template<AgreeSetsGenMethod method = AgreeSetsGenMethod::kUsingVectorOfIDSets>
    SetOfAgreeSets genAgreeSets() const;
    SetOfVectors genPLIMaxRepresentation() const;
    AgreeSet getAgreeSet(int const tuple1_index, int const tuple2_index) const;

private:
    void calculateSupersets(SetOfVectors& max_representation,
                            std::deque<std::vector<int>> const& partition) const;
    ColumnLayoutRelationData const* const relation_;
};
