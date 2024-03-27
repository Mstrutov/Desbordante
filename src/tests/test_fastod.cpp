#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "algorithms/algo_factory.h"
#include "algorithms/od/fastod/fastod.h"
#include "algorithms/od/fastod/hashing/hashing.h"
#include "all_csv_configs.h"
#include "config/names.h"
#include "config/time_limit/type.h"
#include "csv_config_util.h"

namespace tests {

namespace {

size_t RunFastod(CSVConfig const& csv_config) {
    using namespace config::names;

    algos::StdParamsMap params{{kCsvConfig, csv_config}};
    std::unique_ptr<algos::Fastod> fastod = algos::CreateAndLoadAlgorithm<algos::Fastod>(params);

    auto [ods_asc, ods_desc, ods_simple] = fastod->Discover();

    std::vector<algos::fastod::AscCanonicalOD> ods_asc_sorted = ods_asc;
    std::vector<algos::fastod::DescCanonicalOD> ods_desc_sorted = ods_desc;
    std::vector<algos::fastod::SimpleCanonicalOD> ods_simple_sorted = ods_simple;

    std::sort(ods_asc_sorted.begin(), ods_asc_sorted.end());
    std::sort(ods_desc_sorted.begin(), ods_desc_sorted.end());
    std::sort(ods_simple_sorted.begin(), ods_simple_sorted.end());

    size_t ods_asc_sorted_hash = algos::fastod::hashing::CombineHashes(ods_asc_sorted);
    size_t ods_desc_sorted_hash = algos::fastod::hashing::CombineHashes(ods_desc_sorted);
    size_t ods_simple_sorted_hash = algos::fastod::hashing::CombineHashes(ods_simple_sorted);

    std::vector<size_t> od_hashes = {ods_asc_sorted_hash, ods_desc_sorted_hash,
                                     ods_simple_sorted_hash};
    size_t result_hash = algos::fastod::hashing::CombineHashes(od_hashes);

    return result_hash;
}

}  // namespace

class FastodResultHashTest : public ::testing::TestWithParam<CSVConfigHash> {};

TEST_P(FastodResultHashTest, CorrectnessTest) {
    CSVConfigHash csv_config_hash = GetParam();
    size_t actual_hash = RunFastod(csv_config_hash.config);
    EXPECT_EQ(actual_hash, csv_config_hash.hash);
}

INSTANTIATE_TEST_SUITE_P(
        TestFastodSuite, FastodResultHashTest,
        ::testing::Values(CSVConfigHash{kOdTestNormOd, 8741296102670149192ULL},
                          CSVConfigHash{kOdTestNormSmall2x3, 14827049072319306073ULL},
                          CSVConfigHash{kOdTestNormSmall3x3, 66466490561337ULL},
                          CSVConfigHash{kOdTestNormAbalone, 14398696798633970055ULL},
                          CSVConfigHash{kOdTestNormBalanceScale, 11093822414574ULL},
                          CSVConfigHash{kOdTestNormBreastCancerWisconsin, 4334402279000540119ULL},
                          CSVConfigHash{kOdTestNormEchocardiogram, 2243402441338221665ULL},
                          CSVConfigHash{kOdTestNormHorse10c, 1462534374501425106ULL},
                          CSVConfigHash{kOdTestNormIris, 11093822414574ULL},
                          CSVConfigHash{kBernoulliRelation, 6518269127574092257ULL},
                          CSVConfigHash{kTestFD, 15333753345229147120ULL},
                          CSVConfigHash{kWdcAstrology, 723643032648123806ULL},
                          CSVConfigHash{kWdcGame, 3164616462792843131ULL},
                          CSVConfigHash{kWdcPlanetz, 3164616455022529293ULL},
                          CSVConfigHash{kWdcSymbols, 2211268401046792ULL},
                          CSVConfigHash{kNeighbors10k, 11706974185824900569ULL},
                          CSVConfigHash{kNeighbors100k, 11706974185824900569ULL},
                          CSVConfigHash{kAbalone, 13440043079221534278ULL},
                          CSVConfigHash{kIris, 386492228314919716ULL},
                          CSVConfigHash{kBreastCancer, 10457518087798149718ULL}));

}  // namespace tests
