//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#pragma once

#include "Vertical.h"

#include <boost/dynamic_bitset.hpp>

using boost::dynamic_bitset, std::string;

//useless class - think about deprecation
/*class ColumnCombination : public Vertical {

public:
    ColumnCombination(dynamic_bitset<> column_indices_, shared_ptr<RelationalSchema> schema);
    explicit ColumnCombination(Vertical&& vertical) : Vertical(vertical) {}     //or const&??
    string ToString() override ;
};*/
