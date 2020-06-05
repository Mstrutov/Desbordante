//#include "Pyro.h"

#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "algorithms/Pyro.h"

using ::testing::ContainerEq, ::testing::Eq;
using namespace std;

std::string get_selfpath();

//TODO: wrong file => WARNING
TEST(PyroTester, first) {

    //string path = "/home/maxim/Study/Metanome-coding/metanome_git/src/tests";
    auto path = fs::current_path().append("inputData").append("CIPublicHighway50k.csv");
    //CSVParser csvParser(path + "/TestTane.csv");
    //cout << path;
    vector<long> results;
    Pyro algoInstance(path);
    algoInstance.execute();


//    for (auto res : results)
//        cout << res << " ";
    cout << endl;
}