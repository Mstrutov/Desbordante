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
    vector<string> datasets = {"CIPublicHighway50k.csv",
                               //"CIPublicHighway100k.csv",
                               //"iowa50k.csv"//,
                               //"neighbors100k.csv",
                               /*"WDC_astronomical.csv",
                               "WDC_age.csv",
                               "WDC_appearances.csv",
                               "WDC_astrology.csv",
                               "WDC_game.csv",
                               "WDC_science.csv",
                               "WDC_symbols.csv",
                               "WDC_kepler.csv",*/
                               //"adult.csv"
                               };
    // crashes if \r is the newline character
    // string data = "CIPublicHighway50k.csv";
    auto path = fs::current_path().append("inputData");
    //CSVParser csvParser(path + "/TestTane.csv");
    //cout << path;
    try {
        vector<vector<double>> results(1);
        for (int j = 0; j < 1; j++) {
            for (auto & dataset : datasets) {
                Pyro algoInstance(path / dataset, ',', true);
                results[0].push_back(algoInstance.execute());
                std::cout << algoInstance.resultsToJSONString();
            }
        }
        /*Pyro algoInstance(path / "adult.csv", ';', false);
        results.push_back(algoInstance.execute());*/
        cout << "==============PYRO RESULTS============\n";
        for (auto const& algo : results) {
            for (auto time : algo) {
                cout << time << ' ';
            }
        }
        cout << "\n=====================================\n";
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

/*TEST(PyroPerformanceTest, first) {
    vector<tuple<string, char, bool>> datasets = {
            //{"CIPublicHighway50k.csv", ',', true},
            //{"CIPublicHighway100k.csv", ',', true},
            //{"CIPublicHighway200k.csv", ',', true},
            {"adult.csv", ';', false},
            {"breast_cancer.csv", ',', true},
            {"CIPublicHighway.csv", ',', true},
            {"EpicMeds.csv", '|', true},
            {"EpicVitals.csv", '|', true},
            {"iowa1kk.csv", ',', true},
            //{"iowa100k.csv", ',', true},
            {"LegacyPayors.csv", '|', true},
            {"neighbors100k.csv", ',', true}
            //{"digits.csv", ',', true},
    };
    int seed = 1000;
    // crashes if \r is the newline character
    // string data = "CIPublicHighway50k.csv";
    auto path = fs::current_path().append("inputData");
    //CSVParser csvParser(path + "/TestTane.csv");
    //cout << path;
    try {
        int testsNum = 10;
        vector<vector<double>> results(testsNum);
        for (int j = 0; j < testsNum; j++) {
            for (auto& dataset : datasets) {
                cout << "===============DISCOVERING IN " << std::get<0>(dataset) << "===========" << endl;
                Pyro algoInstance(path / std::get<0>(dataset),
                                  std::get<1>(dataset),
                                  std::get<2>(dataset), seed);
                results[j].push_back(algoInstance.execute());
            }
        }
        //Pyro algoInstance(path / "adult.csv", ';', false);
        //results.push_back(algoInstance.execute());
        cout << "==============PYRO RESULTS============\n";
        for (auto const& algo : results) {
            for (auto time : algo) {
                cout << time << ' ';
            }
        }
        cout << "\n=====================================\n";
    }
    catch (std::runtime_error& e) {
        cout << "Exception raised in test: " << e.what() << endl;
        FAIL();
    }

//    for (auto res : results)
//        cout << res << " ";
    cout << endl;
    SUCCEED();
}*/
