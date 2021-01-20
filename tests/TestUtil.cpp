#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "ListAgreeSetSample.h"

TEST(testingBitsetToLonglong, first){
    size_t encoded_num = 1254;
    boost::dynamic_bitset<> simple_bitset{20, encoded_num};

    auto res_vector = *ListAgreeSetSample::bitSetToLongLongVector(simple_bitset);
    ASSERT_EQ(1, res_vector.size());
    for (auto long_long_repr : res_vector)
        ASSERT_EQ(encoded_num, long_long_repr);
}