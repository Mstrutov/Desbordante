//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#pragma once

#include <string>

#include <boost/dynamic_bitset.hpp>
#include <utility>

#include "RelationalSchema.h"

using std::string, boost::dynamic_bitset;

class Column {
friend RelationalSchema;

private:
    string name;
    unsigned int index;
    std::weak_ptr<RelationalSchema> schema;

public:
    Column(std::shared_ptr<RelationalSchema> schema, string name, int index) :
            name(std::move(name)),
            index(index),
            schema(schema) {}
    explicit operator Vertical() const;
    unsigned int getIndex() const { return index; }
    string getName() const { return name; }
    // TODO: const pointer
    std::shared_ptr<RelationalSchema> getSchema() const { return schema.lock(); }
    string toString() const { return "[" + name + "]";}
    std::string toIndicesString() const { return std::to_string(index); }
    explicit operator std::string() const { return toString(); }
    bool operator==(const Column& rhs);
};