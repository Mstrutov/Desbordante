//
// Created by maxim on 18.10.2019.
//

#include "algorithms/TaneX.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>

using ::testing::ContainerEq, ::testing::Eq;
using namespace std;

std::string get_selfpath();

extern string path;

//TODO: wrong file => WARNING
TEST(TaneTester, first){
    Tane algoInstance(path + "/neighbors100k.csv");
    //Tane algoInstance(path + "/CIPublicHighway50k.csv");
    algoInstance.execute();
}
