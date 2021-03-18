//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "Column.h"

class Vertical {
private:
    //Vertical(shared_ptr<RelationalSchema>& relSchema, int indices);

    //TODO: unique_ptr<columnIndices> if this is big
    boost::dynamic_bitset<> columnIndices;
    RelationalSchema const* schema;

public:

    static std::unique_ptr<Vertical> emptyVertical(RelationalSchema const* relSchema);

    Vertical(RelationalSchema const* relSchema, boost::dynamic_bitset<> indices);
    Vertical() = default;

    explicit Vertical(Column const& col);


    Vertical(Vertical const& other) = default;
    Vertical& operator=(const Vertical& rhs) = default;
    Vertical(Vertical&& other) = default;
    Vertical& operator=(Vertical&& rhs) = default;
    ~Vertical() = default;
    bool operator==(Vertical const& other) const { return columnIndices == other.columnIndices; }

    boost::dynamic_bitset<> getColumnIndices() const { return columnIndices; }
    RelationalSchema const* getSchema() const { return schema; }

    bool contains(Vertical const& that) const;
    bool contains(Column const& that) const;
    bool intersects(Vertical const& that) const;
    Vertical Union(Vertical const& that) const;
    Vertical Union(Column const& that) const;
    Vertical project(Vertical const& that) const;
    Vertical without (Vertical const & that) const;
    Vertical without(Column const& that) const;
    Vertical invert() const;
    Vertical invert(Vertical const& scope) const;

    unsigned int getArity() const { return columnIndices.count(); }
    std::vector<Column const*> getColumns() const;
    std::vector<Vertical> getParents() const;

    std::string toIndicesString() const;
    explicit operator std::string() const { return toIndicesString(); }
};
