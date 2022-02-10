#pragma once

#include "Itemset.h"

class AR {
private:
    Itemset left;   //antecedent
    Itemset right;  //consequent
    //double conf;  //нужно ли???
    //TODO может быть не индексы, а сразу строки?
public:
    AR(Itemset const& left, Itemset const& right)
        : left(left), right(right)
    {}

    AR(Itemset && left, Itemset && right)
        :left(std::move(left)), right(std::move(right)) {}
};