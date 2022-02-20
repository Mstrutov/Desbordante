#pragma once

#include "Itemset.h"

class AR {
public:
    std::vector<unsigned> left;   //antecedent
    std::vector<unsigned> right;  //consequent
    double confidence = -1;
    //TODO может быть не индексы, а сразу строки?

    AR() = default;
    AR(std::vector<unsigned> && left, std::vector<unsigned> && right, double confidence)
        :left(std::move(left)), right(std::move(right)), confidence(confidence) {}
};