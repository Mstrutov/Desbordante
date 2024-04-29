#pragma once

#include <string>
#include <utility>

#include <boost/dynamic_bitset.hpp>

#include "column_index.h"

class RelationalSchema;
class Vertical;

class Column {
    friend RelationalSchema;

public:
    using IndexType = model::ColumnIndex;

private:
    std::string name_;
    IndexType index_;
    RelationalSchema const* schema_;

public:
    Column(RelationalSchema const* schema, std::string name, IndexType index)
        : name_(std::move(name)), index_(index), schema_(schema) {}

    IndexType GetIndex() const {
        return index_;
    }

    std::string GetName() const {
        return name_;
    }

    RelationalSchema const* GetSchema() const {
        return schema_;
    }

    std::string ToIndicesString() const {
        return std::to_string(index_);
    }

    std::string ToString() const {
        return "[" + name_ + "]";
    }

    explicit operator std::string() const {
        return ToString();
    }

    explicit operator Vertical() const;

    /* We consider the lhs column to be less than rhs if
     * lhs.index > rhs.index
     */
    bool operator<(Column const& rhs) const;
    bool operator==(Column const& rhs) const;
    bool operator!=(Column const& rhs) const;
    bool operator>(Column const& rhs) const;
};
