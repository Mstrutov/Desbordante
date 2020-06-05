//
// Created by maxim on 18.10.2019.
//

#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "algorithms/TaneX.h"

using ::testing::ContainerEq, ::testing::Eq;
using namespace std;

std::string get_selfpath();

extern string path;

//TODO: wrong file => WARNING
//TODO: NO ABSOLUTE PATHS!!
TEST(TaneTester, first) {

    //string path = "/home/maxim/Study/Metanome-coding/metanome_git/src/tests";
    string path = fs::path(get_selfpath()).parent_path().string() + "/inputData";
    //CSVParser csvParser(path + "/TestTane.csv");
    //cout << path;
    vector<long> results;
    Tane algoInstance(path + "/CIPublicHighway50k.csv");
    //Tane algoInstance(path + "/CIPublicHighway50k.csv");
    results.push_back(algoInstance.execute());
    /*for (int i = 0; i < 10; i++) {
        Tane algoInstance(path + "/neighbors50k.csv");
        //Tane algoInstance(path + "/CIPublicHighway50k.csv");
        results.push_back(algoInstance.execute());
    }
    for (int i = 0; i < 10; i++) {
        Tane algoInstance(path + "/CIPublicHighway100k.csv");
        //Tane algoInstance(path + "/CIPublicHighway50k.csv");
        results.push_back(algoInstance.execute());
    }
    for (int i = 0; i < 10; i++) {
        Tane algoInstance(path + "/LegacyOPMeds.csv");
        //Tane algoInstance(path + "/CIPublicHighway50k.csv");
        results.push_back(algoInstance.execute());
    }
    for (int i = 0; i < 10; i++) {
        Tane algoInstance(path + "/LegacyIPMeds.csv");
        //Tane algoInstance(path + "/CIPublicHighway50k.csv");
        results.push_back(algoInstance.execute());
    }
    for (int i = 0; i < 10; i++) {
        Tane algoInstance(path + "/EpicVitals.csv");
        //Tane algoInstance(path + "/CIPublicHighway50k.csv");
        results.push_back(algoInstance.execute());
    }*/

    for (auto res : results)
        cout << res << " ";
    cout << endl;
}