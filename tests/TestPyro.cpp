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
    vector<string> datasets = {"CIPublicHighway50k.csv", "CIPublicHighway100k.csv", "CIPublicHighway200k.csv", "neighbors100k.csv"};
    // string data = "CIPublicHighway50k.csv";
    auto path = fs::current_path().append("inputData").append(datasets[2]);
    //CSVParser csvParser(path + "/TestTane.csv");
    //cout << path;
    try {
        vector<long> results;
        Pyro algoInstance(path);
        algoInstance.execute();
    }
    catch (std::runtime_error& e) {
        cout << "Exception raised in test: " << e.what() << endl;
        FAIL();
    }

//    for (auto res : results)
//        cout << res << " ";
    cout << endl;
    SUCCEED();
}