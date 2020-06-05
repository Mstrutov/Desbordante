//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "ColumnLayoutRelationData.h"

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

    string path = fs::path(get_selfpath()).parent_path().string() + "/inputData/Test1.csv";
    CSVParser csvParser(path);
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

    string path = fs::path(get_selfpath()).parent_path().string()  + "/inputData/Test1.csv";
    CSVParser csvParser(path);
    auto test = ColumnLayoutRelationData::createFrom(csvParser, false);
    auto columnData = test->getColumnData(0);
    auto index = columnData->getPositionListIndex()->getIndex();
    ASSERT_THAT(ans, ContainerEq(index));
}

TEST(pliIntersectChecker, first){
    deque<vector<int>> ans = {
            {2, 5}
    };

    string path = fs::path(get_selfpath()).parent_path().string() + "/inputData";
    CSVParser csvParser1(path + "/ProbeTest1.csv");
    CSVParser csvParser2(path + "/ProbeTest2.csv");
    auto test1 = ColumnLayoutRelationData::createFrom(csvParser1, false);
    auto test2 = ColumnLayoutRelationData::createFrom(csvParser2, false);

    auto pli1 = test1->getColumnData(0)->getPositionListIndex();
    auto pli2 = test2->getColumnData(0)->getPositionListIndex();

    auto intersection = pli1->intersect(pli2);
    int i = 0;
    //auto index = columnData->getPositionListIndex()->getIndex();
    //ASSERT_THAT(intersection->getIndex(), ContainerEq(ans));
}