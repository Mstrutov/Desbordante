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
    kUsingVectorOfIDSets = 0, /*< Generates agree sets using identifier sets.
                               *  Algorithm works as follows:
                               *  1. Generates maximal representation
                               *     (check out the `kUsingMCAndGetAgreeSet` description).
                               *  2. Fills vector<IdentifierSet> with identifier sets by
                               *     iterating over each cluster in max representation.
                               *     In order to avoid adding identifier set of the same tuple
                               *     twice (the same tuple can be in different clusters), the
                               *     set of ids of the already added tuples is used.
                               *  3. Iterates over all pairs of identifier sets from vector.
                               *  4. Gets agree set for current pair by intersecting.
                               */
    kUsingMapOfIDSets       , /*< Metanome approach.
                               *  Generates agree sets using identifier sets.
                               *  Algorithm works as follows:
                               *  1. Generates maximal representation
                               *     (check out the `kUsingMCAndGetAgreeSet` description).
                               *  2. Fills map<int, IdentifierSet> with <tuple index, tuple idset>
                               *     pairs by iterating over each cluster in max representation.
                               *  3. Iterates over all pairs of tuples from each cluster
                               *     of maximal representation.
                               *  4. Gets agree set for current pair of tuples by intersecting
                               *     their identifier sets.
                               */
    kUsingGetAgreeSet       , /*< The most naive (so the slowest) way to generate agree sets.
                               *  Generates agree set for all pairs of tuples that
                               *  can form agree set. Tuples can form agree set if
                               *  they have the same value in at least one attribute.
                               *  Such tuples are in the same cluster of pli.
                               *  So algorithm works as follows:
                               *  1. Iterates over all table attributes.
                               *  2. Iterates over all pairs of tuples from each cluster
                               *     of current attribute pli.
                               *  3. Gets agree set for current pair using getAgreeSet.
                               */
    kUsingMCAndGetAgreeSet  , /*< In `kUsingGetAgreeSet` method, the same pair of tuples
                               *  can be processed (passed to getAgreeSet) multiple times.
                               *  This method similar to `kUsingGetAgreeSet`,
                               *  but avoids the above problem.
                               *  It reduces the number of redundant tuple comparisons by
                               *  transforming partitions into its maximal representation.
                               *  And processes tuple pairs from maximal representation clusters.
                               *  Algorithm of maximal representation generation was taken from
                               *  Metanome. For more information about maximal representation
                               *  check out http://www.vldb.org/pvldb/vol8/p1082-papenbrock.pdf
                               */
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
