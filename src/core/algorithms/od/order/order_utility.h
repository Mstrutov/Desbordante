#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/container_hash/hash.hpp>

#include "model/table/column_index.h"
#include "model/table/tuple_index.h"
#include "model/table/typed_column_data.h"

namespace algos::order {
using Node = std::vector<model::ColumnIndex>;
using AttributeList = std::vector<model::ColumnIndex>;
using CandidatePairs = std::vector<std::pair<AttributeList, AttributeList>>;
using Hash = boost::hash<AttributeList>;
using CandidateSets =
        std::unordered_map<AttributeList, std::unordered_set<AttributeList, Hash>, Hash>;
using OrderDependencies =
        std::unordered_map<AttributeList, std::unordered_set<AttributeList, Hash>, Hash>;

struct IndexedByteData {
    model::TupleIndex index;
    std::byte const* data;
};

void PrintOD(AttributeList const& lhs, AttributeList const& rhs);
std::vector<AttributeList> GetPrefixes(Node const& node);
AttributeList MaxPrefix(AttributeList const& attribute_list);
bool InUnorderedMap(OrderDependencies const& map, AttributeList const& lhs,
                    AttributeList const& rhs);
bool AreDisjoint(AttributeList const& a, AttributeList const& b);
bool StartsWith(AttributeList const& rhs_candidate, AttributeList const& rhs);
std::unordered_set<model::TupleIndex> GetNullIndices(
        std::vector<model::TypedColumnData> const& data);
std::vector<IndexedByteData> GetIndexedByteData(
        model::TypedColumnData const& data, std::unordered_set<model::TupleIndex> const& null_rows);
}  // namespace algos::order
