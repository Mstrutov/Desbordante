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
TEST(TaneTester, first) {
    vector<string> datasets = {"CIPublicHighway50k.csv",
                               "CIPublicHighway100k.csv",
                               "CIPublicHighway200k.csv",
                               "neighbors100k.csv",
                               "WDC_astronomical.csv",
                               "WDC_age.csv",
                               "WDC_appearances.csv",
                               "WDC_astrology.csv",
                               "WDC_game.csv",
                               "WDC_science.csv",
                               "WDC_symbols.csv",
                               "WDC_kepler.csv"};
    auto path = fs::current_path().append("inputData");
    vector<long> results;
    try {
        /*for (int i = 0; i < datasets.size(); i++) {
            Tane algoInstance(path / datasets[i]);
            results.push_back(algoInstance.execute());
        }*/
        cout << "==============TANE RESULTS============\n";
        for (auto time : results) {
            cout << time << ' ';
        }
        cout << "\n=====================================\n";
    }
    catch (std::runtime_error& e) {
        cout << "Exception raised in test: " << e.what() << endl;
        FAIL();
    }

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
    SUCCEED();
}