#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "Pyro.h"
#include "algorithms/TaneX.h"

using ::testing::ContainerEq, ::testing::Eq;

using std::string, std::vector;

TEST(TaneTester, againstPyro) {
    vector<string> datasets = {"CIPublicHighway50k.csv",
                               "neighbors100k.csv",
                               "WDC_astronomical.csv",
                               "WDC_age.csv",
                               "WDC_appearances.csv",
                               "WDC_astrology.csv",
                               "WDC_game.csv",
                               "WDC_science.csv",
                               "WDC_symbols.csv",
                               "WDC_kepler.csv"};
    vector<char> separators = {',', ',', ',', ',', ',', ',', ',', ',', ',', ','};
    vector<bool> header_presence = {true, true, true, true, true, true, true, true, true, true};
    auto path = fs::current_path().append("inputData");

    try {
        for (size_t i = 0; i < datasets.size(); i++) {
            auto tane = Tane(path / datasets[i], separators[i], header_presence[i]);
            tane.execute();
            std::string resultsTane = tane.getJsonFDs();
            auto pyro = Pyro(path / datasets[i], separators[i], header_presence[i], 0, 0, -1);
            pyro.execute();
            std::string resultsPyro = pyro.FDAlgorithm::getJsonFDs();
            ASSERT_EQ(resultsPyro, resultsTane) << "Tane and Pyro yield different results at " << datasets[i];
        }
    }
    catch (std::runtime_error& e) {
        std::cout << "Exception raised in test: " << e.what() << std::endl;
        FAIL();
    }
    SUCCEED();
}
