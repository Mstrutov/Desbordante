#include <filesystem>
#include <iostream>
#include <string>

#include <boost/container_hash/hash.hpp>
#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/od/order/order.h"
#include "all_csv_configs.h"
#include "config/names.h"
#include "csv_config_util.h"

namespace tests {

namespace fs = std::filesystem;

using OD = algos::order::OrderDependencies;

class OrderTest : public ::testing::Test {
public:
    static std::unique_ptr<algos::order::Order> CreateOrderInstance(CSVConfig const& info) {
        using namespace config::names;
        return algos::CreateAndLoadAlgorithm<algos::order::Order>({{kCsvConfig, info}});
    }
};

TEST_F(OrderTest, SmallDataset) {
    auto a = CreateOrderInstance(kODnorm6);
    a->Execute();
    OD actual = a->GetValidODs();

    OD expected;
    expected[{0}] = {{1}, {3}, {4}};
    expected[{1}] = {{4}, {0, 2}, {0, 5}, {3, 2}, {3, 5}};
    expected[{2}] = {{5}};
    expected[{3}] = {{4}, {0, 5}, {1, 5}};
    expected[{2, 0}] = {{5}};
    expected[{2, 1}] = {{5}};
    expected[{2, 3}] = {{5}};
    expected[{2, 4}] = {{5}};
    expected[{4}] = {{0, 5}, {1, 3}, {1, 5}, {3, 1}, {3, 2}, {3, 5}, {1, 2, 3}, {0, 2, 3}};
    expected[{0, 2}] = {{1}, {4}, {3, 1}, {3, 5}};
    expected[{0, 5}] = {{4}, {1, 3}, {3, 1}, {3, 2}, {1, 2, 3}};
    expected[{1, 2}] = {{4}, {0, 5}, {3, 5}};
    expected[{1, 3}] = {{4}, {0, 5}};
    expected[{1, 5}] = {{4}, {3, 2}, {0, 2, 3}};
    expected[{3, 1}] = {{4}, {0, 5}};
    expected[{3, 2}] = {{4}, {0, 5}, {1, 5}};
    expected[{3, 5}] = {{4}};
    expected[{5}] = {{2, 3}, {2, 4}, {2, 1, 3}};
    expected[{1, 2, 3}] = {{4}, {0, 5}};
    expected[{0, 2, 3}] = {{4}, {1, 5}};
    expected[{2, 1, 3}] = {{5}};

    EXPECT_EQ(expected, actual);
}

TEST_F(OrderTest, BigWithDifferentTypes) {
    auto a = CreateOrderInstance(kNeighbors10k);
    a->Execute();
    OD actual = a->GetValidODs();

    OD expected;
    expected[{5}] = {{0}, {1}, {2}, {3}, {4}, {6}};

    EXPECT_EQ(expected, actual);
}

TEST_F(OrderTest, OrderTest_Epic_Test) {
    auto a = CreateOrderInstance(kEpicMeds);
    a->Execute();
    OD actual = a->GetValidODs();

    OD expected;
    expected[{9}] = {{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}};

    EXPECT_EQ(expected, actual);
}

}  // namespace tests
