//#include "Pyro.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iostream>
#include "algorithms/Pyro.h"

using ::testing::ContainerEq, ::testing::Eq;
using namespace std;

std::string get_selfpath();

//TODO: wrong file => WARNING
//TODO: NO ABSOLUTE PATHS!!
TEST(PyroTester, first) {

    //string path = "/home/maxim/Study/Metanome-coding/metanome_git/src/tests";
    string path = fs::path(get_selfpath()).parent_path().string() + "/inputData";
    //CSVParser csvParser(path + "/TestTane.csv");
    //cout << path;
    vector<long> results;
    Pyro algoInstance(path + "/CIPublicHighway50k.csv");
    algoInstance.execute();


//    for (auto res : results)
//        cout << res << " ";
    cout << endl;
}