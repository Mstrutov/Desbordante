#pragma once

#include <vector>
#include <memory>

#include "ColumnData.h"
#include "Vertical.h"
#include "ColumnLayoutRelationData.h"

/* TODO: consider using Vertical + vector<int> indices
 * instead of vector<IdentitfierSetValue>. Probably it will be faster
 * but also maybe will hurt readability.
 * UPD: Did not notice any difference in perfmormance
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

    ColumnLayoutRelationData const* relation_;
    std::vector<IdentifierSetValue> data_;
    int tuple_index_;
};