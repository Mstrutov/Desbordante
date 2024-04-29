//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/optional.hpp>

#include "model/table/column.h"
#include "model/table/idataset_stream.h"
#include "model/table/vertical.h"
#include "util/bitset_utils.h"

class RelationalSchema {
private:
    std::vector<Column> columns_;
    std::string name_;

public:
    std::unique_ptr<Vertical> empty_vertical_;

    RelationalSchema(std::string name, std::vector<std::string> column_names);

    RelationalSchema(RelationalSchema const& other) = delete;
    RelationalSchema& operator=(RelationalSchema const& rhs) = delete;
    RelationalSchema(RelationalSchema&& other) noexcept = default;
    RelationalSchema& operator=(RelationalSchema&& rhs) noexcept = default;

    static std::unique_ptr<RelationalSchema> CreateFrom(model::IDatasetStream& table);

    std::string GetName() const {
        return name_;
    }

    std::vector<Column> const& GetColumns() const {
        return columns_;
    };

    Column const& GetColumn(std::string const& col_name) const;
    Column const& GetColumn(size_t index) const;
    size_t GetNumColumns() const;
    Vertical GetVertical(boost::dynamic_bitset<> indices) const;

    template <typename Container>
    boost::dynamic_bitset<> IndicesToBitset(Container const& indices) const;
    template <typename ForwardIt>
    boost::dynamic_bitset<> IndicesToBitset(ForwardIt begin, ForwardIt end) const;

    std::unordered_set<Vertical> CalculateHittingSet(
            std::vector<Vertical> verticals,
            boost::optional<std::function<bool(Vertical const&)>> pruning_function) const;

    ~RelationalSchema();

    friend inline bool operator==(RelationalSchema const& l, RelationalSchema const& r);
};

inline bool operator==(RelationalSchema const& l, RelationalSchema const& r) {
    return (l.name_ == r.name_ && l.columns_.size() == r.columns_.size());
}

inline bool operator!=(RelationalSchema const& l, RelationalSchema const& r) {
    return !(l == r);
}

template <typename Container>
boost::dynamic_bitset<> RelationalSchema::IndicesToBitset(Container const& indices) const {
    return util::IndicesToBitset(indices, GetNumColumns());
}

template <typename ForwardIt>
boost::dynamic_bitset<> RelationalSchema::IndicesToBitset(ForwardIt begin, ForwardIt end) const {
    return util::IndicesToBitset(begin, end, GetNumColumns());
}
