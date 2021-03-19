#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "algorithms/Pyro.h"
#include "Datasets.h"

using ::testing::ContainerEq, ::testing::Eq;
using std::vector, std::string, std::cout, std::endl;



class PyroTest : public LightDatasets, public HeavyDatasets, public ::testing::Test {};

TEST_F(PyroTest, ReturnsSameHashOnHeavyDatasets) {

    auto path = std::filesystem::current_path().append("inputData");

    try {
        for (auto dataset : HeavyDatasets::datasets) {
            auto pyro = Pyro(path / dataset.name, dataset.separator, dataset.header_presence, 0, 0, -1);
            pyro.execute();
            EXPECT_EQ(pyro.FDAlgorithm::fletcher16(), dataset.hash)
                << "Pyro result hash changed for " << dataset.name;
        }
    }
    catch (std::runtime_error& e) {
        std::cout << "Exception raised in test: " << e.what() << std::endl;
        FAIL();
    }
    SUCCEED();
}
