#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "ColumnLayoutRelationData.h"
#include "ListAgreeSetSample.h"

using ::testing::ContainerEq, ::testing::Eq;
using namespace std;

std::string get_selfpath();

TEST(pliChecker, first){
    deque<vector<int>> ans = {
            {0, 2, 8, 11},
            {1, 5, 9},
            {4, 14},
            {6, 7, 18},
            {10, 17} //null
    };

    auto path = fs::current_path().append("inputData").append("Test1.csv");
    deque<vector<int>> index;
    try {
        CSVParser csvParser(path);
        auto test = ColumnLayoutRelationData::createFrom(csvParser, true);
        auto columnData = test->getColumnData(0);
        index = columnData->getPositionListIndex()->getIndex();
    }
    catch (std::runtime_error& e) {
        cout << "Excepion raised in test: " << e.what() << endl;
        FAIL();
    }
    ASSERT_THAT(ans, ContainerEq(index));
}

TEST(pliChecker, second){
    deque<vector<int>> ans = {
            {0, 2, 8, 11},
            {1, 5, 9},
            {4, 14},
            {6, 7, 18},
    };
    deque<vector<int>> index;
    try {
        auto path = fs::current_path().append("inputData").append("Test1.csv");
        CSVParser csvParser(path);
        auto test = ColumnLayoutRelationData::createFrom(csvParser, false);
        auto columnData = test->getColumnData(0);
        index = columnData->getPositionListIndex()->getIndex();
    }
    catch (std::runtime_error& e) {
        cout << "Excepion raised in test: " << e.what() << endl;
        FAIL();
    }
    ASSERT_THAT(ans, ContainerEq(index));
}

TEST(pliIntersectChecker, first){
    deque<vector<int>> ans = {
            {2, 5}
    };
    shared_ptr<PositionListIndex> intersection;

    try {
        auto path = fs::current_path().append("inputData");
        CSVParser csvParser1(path / "ProbeTest1.csv");
        CSVParser csvParser2(path / "ProbeTest2.csv");

        auto test1 = ColumnLayoutRelationData::createFrom(csvParser1, false);
        auto test2 = ColumnLayoutRelationData::createFrom(csvParser2, false);
        auto pli1 = test1->getColumnData(0)->getPositionListIndex();
        auto pli2 = test2->getColumnData(0)->getPositionListIndex();

        intersection = pli1->intersect(pli2);
    }
    catch (std::runtime_error& e) {
        cout << "Excepion raised in test: " << e.what() << endl;
        FAIL();
    }
    ASSERT_THAT(intersection->getIndex(), ContainerEq(ans));
}

TEST(testingBitsetToLonglong, first){
    size_t encoded_num = 1254;
    boost::dynamic_bitset<> simple_bitset{20, encoded_num};

    auto res_vector = *ListAgreeSetSample::bitSetToLongLongVector(simple_bitset);
    ASSERT_EQ(1, res_vector.size());
    for (auto long_long_repr : res_vector)
        ASSERT_EQ(encoded_num, long_long_repr);
}