#include <filesystem>

#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/fd/fastfds/fastfds.h"
#include "algorithms/fd/fd_algorithm.h"
#include "algorithms/fd/tane/tane.h"
#include "config/names.h"
#include "datasets.h"
#include "model/table/relational_schema.h"

using algos::FastFDs, algos::Tane;
using config::InputTable;

InputTable GetInputTable(std::filesystem::path const& path, char separator = ',',
                         bool has_header = true) {
    return std::make_unique<CSVParser>(path, separator, has_header);
}

template <typename T>
class MaxLHSOptionTest : public ::testing::Test {
protected:
    static std::unique_ptr<algos::FDAlgorithm> CreateAlgoInstance(unsigned int max_lhs,
                                                                  bool EqNullsType,
                                                                  std::filesystem::path const& path,
                                                                  char separator = ',',
                                                                  bool hasHeader = true) {
        using namespace config::names;

        algos::StdParamsMap params{
                {kCsvPath, path},        {kSeparator, separator},
                {kHasHeader, hasHeader}, {kTable, GetInputTable(path, separator, hasHeader)},
                {kMaximumLhs, max_lhs},  {kEqualNulls, EqNullsType}};

        return algos::CreateAndLoadAlgorithm<T>(params);
    }
};

TYPED_TEST_SUITE_P(MaxLHSOptionTest);

TYPED_TEST_P(MaxLHSOptionTest, MaxLHSOptionWork) {
    unsigned int max_lhs = 1;
    auto const path = test_data_dir / "TestFD.csv";
    auto algorithm = TestFixture::CreateAlgoInstance(max_lhs, true, path);
    algorithm->Execute();
    std::list<FD> result_fds_list = algorithm->FdList();
    for (auto& fd : result_fds_list) {
        ASSERT_TRUE(fd.GetLhs().GetArity() <= max_lhs);
    }
}

REGISTER_TYPED_TEST_SUITE_P(MaxLHSOptionTest, MaxLHSOptionWork);

using Algorithms = ::testing::Types<algos::Tane, algos::FastFDs>;

INSTANTIATE_TYPED_TEST_SUITE_P(MaxLHSOptionTest, MaxLHSOptionTest, Algorithms);