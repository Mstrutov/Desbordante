
#include <iostream>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "algorithms/Pyro.h"
#include "Datasets.h"

using ::testing::ContainerEq, ::testing::Eq;
using std::vector, std::string, std::cout, std::endl;



class PyroTest : public LightDatasets, public HeavyDatasets, public ::testing::Test {};

TEST_F(PyroTest, ReturnsSameHashOnHeavyDatasets) {

    auto path = fs::current_path().append("inputData");

    try {
        for (size_t i = 0; i < HeavyDatasets::datasetQuantity(); i++) {
            auto pyro = Pyro(path / HeavyDatasets::dataset(i), HeavyDatasets::separator(i),
                             HeavyDatasets::hasHeader(i), 0, 0, -1);
            pyro.execute();
            EXPECT_EQ(pyro.FDAlgorithm::fletcher16(), HeavyDatasets::hash(i))
                << "Pyro result hash changed for " << HeavyDatasets::dataset(i);
        }
    }
    catch (std::runtime_error& e) {
        std::cout << "Exception raised in test: " << e.what() << std::endl;
        FAIL();
    }
    SUCCEED();
}
