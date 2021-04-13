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

using boost::dynamic_bitset, std::string, std::weak_ptr, std::shared_ptr, std::vector;

class Vertical {
private:
    //Vertical(shared_ptr<RelationalSchema>& relSchema, int indices);

    dynamic_bitset<> columnIndices;
    weak_ptr<RelationalSchema> schema;

public:
    Vertical(shared_ptr<RelationalSchema> relSchema, dynamic_bitset<> const & indices);
    Vertical() = default;
    explicit Vertical(Column & col);

    //Vertical(Vertical& other) = default;
    Vertical(Vertical const& other) = default;
    Vertical& operator=(const Vertical& rhs) = default;
    Vertical(Vertical&& other) = default;
    Vertical& operator=(Vertical&& rhs) = default;
    virtual ~Vertical() = default;
    bool operator==(Vertical const& other) const { return columnIndices == other.columnIndices; }

    dynamic_bitset<> getColumnIndices() const { return columnIndices; }
    shared_ptr<RelationalSchema> getSchema() { return schema.lock(); }
    bool contains(Vertical const& that) const;
    bool intersects(Vertical const& that) const;
    std::shared_ptr<Vertical> Union(Vertical const& that) const;
    std::shared_ptr<Vertical> project(Vertical const& that) const;
    std::shared_ptr<Vertical> without (Vertical const & that) const;
    std::shared_ptr<Vertical> invert() const;
    std::shared_ptr<Vertical> invert(Vertical const& scope) const;

    static Vertical emptyVertical(shared_ptr<RelationalSchema> relSchema);
    unsigned int getArity() const { return columnIndices.count(); }
    vector<shared_ptr<Column>> getColumns() const;
    vector<shared_ptr<Vertical>> getParents();

    string toString() const;

    // returns a Vertical as a string "[index_1,index_2,...,index_n]"
    std::string toIndicesString() const;
    explicit operator std::string() const { return toString(); }
};
