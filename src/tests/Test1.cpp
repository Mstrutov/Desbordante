//
// Created by kek on 07.08.2019.
//

#include "model/ColumnLayoutRelationData.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>

using ::testing::ContainerEq, ::testing::Eq;
using namespace std;

std::string get_selfpath();

string path = fs::path(get_selfpath()).parent_path().string() + "/inputData";

TEST(pliChecker, first){
    deque<vector<int>> ans = {
            {0, 2, 8, 11},
            {1, 5, 9},
            {4, 14},
            {6, 7, 18},
            {10, 17} //null
    };

    cout << path << endl;
    CSVParser csvParser(path + "/Test1.csv");
    auto test = ColumnLayoutRelationData::createFrom(csvParser, true);
    auto columnData = test->getColumnData(0);
    auto index = columnData->getPositionListIndex()->getIndex();
    ASSERT_THAT(ans, ContainerEq(index));
}

TEST(pliChecker, second){
    deque<vector<int>> ans = {
            {0, 2, 8, 11},
            {1, 5, 9},
            {4, 14},
            {6, 7, 18},
    };

    CSVParser csvParser(path + "/Test1.csv");
    auto test = ColumnLayoutRelationData::createFrom(csvParser, false);
    auto columnData = test->getColumnData(0);
    auto index = columnData->getPositionListIndex()->getIndex();
    ASSERT_THAT(ans, ContainerEq(index));
}

/*TEST(pliIntersectChecker, first){
    vector<vector<int>> ans = {
            {2, 5}
    };

    CSVParser csvParser1(path + "/ProbeTest1.csv");
    CSVParser csvParser2(path + "/ProbeTest2.csv");
    auto test1 = ColumnLayoutRelationData::createFrom(csvParser1, false);
    auto test2 = ColumnLayoutRelationData::createFrom(csvParser2, false);

    auto pli1 = test1->getColumnData(0)->getPositionListIndex();
    auto pli2 = test2->getColumnData(0)->getPositionListIndex();

    auto intersection = pli1->intersect(pli2);
    int i = 0;
    auto index = columnData->getPositionListIndex()->getIndex();
    ASSERT_THAT(intersection->getIndex(), ContainerEq(ans));
}*/