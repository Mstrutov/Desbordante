//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "column.h"

class RelationalSchema;

class Vertical {
private:
    // Vertical(shared_ptr<RelationalSchema>& relSchema, int indices);

    // TODO: unique_ptr<column_indices_> if this is big
    boost::dynamic_bitset<> column_indices_;
    RelationalSchema const* schema_;

public:
    static std::unique_ptr<Vertical> EmptyVertical(RelationalSchema const* rel_schema);

    Vertical(RelationalSchema const* rel_schema, boost::dynamic_bitset<> indices);
    Vertical() = default;

    explicit Vertical(Column const& col);

    Vertical(Vertical const& other) = default;
    Vertical& operator=(Vertical const& rhs) = default;
    Vertical(Vertical&& other) = default;
    Vertical& operator=(Vertical&& rhs) = default;

    virtual ~Vertical() = default;

    /* @return Returns true if lhs.column_indices_ lexicographically less than
     * rhs.column_indices_ treating bitsets big endian.
     * @brief We do not use directly boost::dynamic_bitset<> operator< because
     * it treats bitsets little endian during comparison and this is not
     * suitable for this case, check out operator< for Columns.
     */
    bool operator<(Vertical const& rhs) const;

    bool operator==(Vertical const& other) const {
        return column_indices_ == other.column_indices_;
    }

    bool operator!=(Vertical const& other) const {
        return column_indices_ != other.column_indices_;
    }

    bool operator>(Vertical const& rhs) const {
        return !(*this < rhs && *this == rhs);
    }

    boost::dynamic_bitset<> GetColumnIndices() const {
        return column_indices_;
    }

    boost::dynamic_bitset<> const& GetColumnIndicesRef() const {
        return column_indices_;
    }

    RelationalSchema const* GetSchema() const {
        return schema_;
    }

    bool Contains(Vertical const& that) const;
    bool Contains(Column const& that) const;
    bool Intersects(Vertical const& that) const;
    Vertical Union(Vertical const& that) const;
    Vertical Union(Column const& that) const;
    Vertical Project(Vertical const& that) const;
    Vertical Without(Vertical const& that) const;
    Vertical Without(Column const& that) const;
    Vertical Invert() const;
    Vertical Invert(Vertical const& scope) const;

    unsigned int GetArity() const {
        return column_indices_.count();
    }

    std::vector<Column const*> GetColumns() const;
    std::vector<unsigned> GetColumnIndicesAsVector() const;
    std::vector<Vertical> GetParents() const;

    std::string ToString() const;
    std::string ToIndicesString() const;

    explicit operator std::string() const {
        return ToString();
    }
};
