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

//TODO: wrong file => WARNING
//TODO: NO ABSOLUTE PATHS!!
TEST(TaneTester, first){

    string path = "/home/maxim/Study/Metanome-coding/metanome_git/src/tests";
    //string path = fs::path(get_selfpath()).parent_path().string();
    //CSVParser csvParser(path + "/TestTane.csv");
    Tane algoInstance(path + "/neighbors100k.csv");
    //Tane algoInstance(path + "/CIPublicHighway50k.csv");
    algoInstance.execute();
}
