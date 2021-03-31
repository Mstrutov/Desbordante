
#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "algorithms/Pyro.h"

using ::testing::ContainerEq, ::testing::Eq;
using std::vector, std::string, std::cout, std::endl;

TEST(PyroTester, heavyDatasetsHash) {
    vector<unsigned int> hashes = {20873, 58641, 32696, 19617, 24082, 55285, 34874, 62210};
    vector<string> datasets = {"adult.csv", "breast_cancer.csv", "CIPublicHighway.csv", "EpicMeds.csv",
                               "EpicVitals.csv", "iowa1kk.csv", "LegacyPayors.csv", "neighbors100k.csv"};
    vector<char> separators = {';', ',', ',', '|', '|', ',', '|', ','};
    vector<bool> hasHeader = {false, true, true, true, true, true, true, true};

    auto path = fs::current_path().append("inputData");

    try {
        for (size_t i = 0; i < datasets.size(); i++) {
            auto pyro = Pyro(path / datasets[i], separators[i], hasHeader[i], 0, 0, -1);
            pyro.execute();
            ASSERT_EQ(pyro.FDAlgorithm::fletcher16(), hashes[i]) << "Pyro result hash changed at " << datasets[i];
        }
    }
    catch (std::runtime_error& e) {
        std::cout << "Exception raised in test: " << e.what() << std::endl;
        FAIL();
    }
    SUCCEED();
}
