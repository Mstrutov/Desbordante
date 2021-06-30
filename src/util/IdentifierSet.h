#pragma once

#include <vector>
#include <memory>

#include "ColumnData.h"
#include "Vertical.h"
#include "ColumnLayoutRelationData.h"

/* Class which represents the relationship between a tuple and
 * all partitions containing it. Given the tuple t, IdentifierSet
 * stores vector of <attribute, cluster_index> pairs, where `cluster_index`
 * is the index of cluster in `attribute` pli the t belongs to.
 * Intersection of two identifier sets is the agree set for appropriate tuples.
 * For more information check out http://www.vldb.org/pvldb/vol8/p1082-papenbrock.pdf
 * TODO: consider using Vertical + vector<int> indices
 * instead of vector<IdentitfierSetValue>. Probably it will be faster
 * but also maybe will hurt readability.
 * UPD: Did not notice any difference in perfmormance.
 */
class IdentifierSet {
public:
    IdentifierSet(IdentifierSet const&) = default;
    IdentifierSet(IdentifierSet&&) = default;
    IdentifierSet(ColumnLayoutRelationData const* relation, int index);

    std::string toString() const;

    // Returns an intersection (agree_set(tuple, other.tuple)) of two IndetifierSets
    Vertical intersect(IdentifierSet const& other) const;
private:
    struct IdentifierSetValue {
        Column const* attribute;
        int cluster_index;
    };

    ColumnLayoutRelationData const* const relation_;
    std::vector<IdentifierSetValue> data_;
    int const tuple_index_;
};
