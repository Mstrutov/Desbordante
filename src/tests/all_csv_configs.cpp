#include "all_csv_configs.h"

#include <string_view>

#include "csv_config_util.h"

namespace tests {

namespace {
/// create `CSVConfig` using relative path to the directory with test data
CSVConfig CreateCsvConfig(std::string_view filename, char separator, bool has_header) {
    return {kTestDataDir / filename, separator, has_header};
}
}  // namespace

CSVConfig const kWdcAstronomical = CreateCsvConfig("WDC_astronomical.csv", ',', true);
CSVConfig const kWdcSymbols = CreateCsvConfig("WDC_symbols.csv", ',', true);
CSVConfig const kWdcScience = CreateCsvConfig("WDC_science.csv", ',', true);
CSVConfig const kWdcSatellites = CreateCsvConfig("WDC_satellites.csv", ',', true);
CSVConfig const kWdcAppearances = CreateCsvConfig("WDC_appearances.csv", ',', true);
CSVConfig const kWdcAstrology = CreateCsvConfig("WDC_astrology.csv", ',', true);
CSVConfig const kWdcGame = CreateCsvConfig("WDC_game.csv", ',', true);
CSVConfig const kWdcKepler = CreateCsvConfig("WDC_kepler.csv", ',', true);
CSVConfig const kWdcPlanetz = CreateCsvConfig("WDC_planetz.csv", ',', true);
CSVConfig const kWdcAge = CreateCsvConfig("WDC_age.csv", ',', true);
CSVConfig const kTestWide = CreateCsvConfig("TestWide.csv", ',', true);
CSVConfig const kAbalone = CreateCsvConfig("abalone.csv", ',', false);
CSVConfig const kIris = CreateCsvConfig("iris.csv", ',', false);
CSVConfig const kAdult = CreateCsvConfig("adult.csv", ';', false);
CSVConfig const kBreastCancer = CreateCsvConfig("breast_cancer.csv", ',', true);
CSVConfig const kCIPublicHighway10k = CreateCsvConfig("CIPublicHighway10k.csv", ',', true);
CSVConfig const kNeighbors10k = CreateCsvConfig("neighbors10k.csv", ',', true);
CSVConfig const kNeighbors50k = CreateCsvConfig("neighbors50k.csv", ',', true);
CSVConfig const kNeighbors100k = CreateCsvConfig("neighbors100k.csv", ',', true);
CSVConfig const kCIPublicHighway700 = CreateCsvConfig("CIPublicHighway700.csv", ',', true);
CSVConfig const kEpicVitals = CreateCsvConfig("EpicVitals.csv", '|', true);
CSVConfig const kEpicMeds = CreateCsvConfig("EpicMeds.csv", '|', true);
CSVConfig const kIowa1kk = CreateCsvConfig("iowa1kk.csv", ',', true);
CSVConfig const kFdReduced30 = CreateCsvConfig("fd-reduced-30.csv", ',', true);
CSVConfig const kFlight1k = CreateCsvConfig("flight_1k.csv", ';', true);
CSVConfig const kPlista1k = CreateCsvConfig("plista_1k.csv", ';', false);
CSVConfig const kLetter = CreateCsvConfig("letter.csv", ',', false);
CSVConfig const kCIPublicHighway = CreateCsvConfig("CIPublicHighway.csv", ',', true);
CSVConfig const kLegacyPayors = CreateCsvConfig("LegacyPayors.csv", '|', true);
CSVConfig const kTestEmpty = CreateCsvConfig("TestEmpty.csv", ',', true);
CSVConfig const kTestSingleColumn = CreateCsvConfig("TestSingleColumn.csv", ',', true);
CSVConfig const kTestLong = CreateCsvConfig("TestLong.csv", ',', true);
CSVConfig const kTestFD = CreateCsvConfig("TestFD.csv", ',', true);
CSVConfig const kTestND = CreateCsvConfig("TestND.csv", ',', true);
CSVConfig const kOdTestNormOd = CreateCsvConfig("od_norm_data/OD_norm.csv", ',', true);
CSVConfig const kOdTestNormSmall2x3 = CreateCsvConfig("od_norm_data/small_2x3.csv", ',', true);
CSVConfig const kOdTestNormSmall3x3 = CreateCsvConfig("od_norm_data/small_3x3.csv", ',', true);
CSVConfig const kOdTestNormAbalone =
        CreateCsvConfig("od_norm_data/metanome/abalone_norm.csv", ',', true);
CSVConfig const kOdTestNormBalanceScale =
        CreateCsvConfig("od_norm_data/metanome/balance-scale_norm.csv", ',', true);
CSVConfig const kOdTestNormBreastCancerWisconsin =
        CreateCsvConfig("od_norm_data/metanome/breast-cancer-wisconsin.csv", ',', true);
CSVConfig const kOdTestNormEchocardiogram =
        CreateCsvConfig("od_norm_data/metanome/echocardiogram_norm.csv", ',', true);
CSVConfig const kOdTestNormHorse10c =
        CreateCsvConfig("od_norm_data/metanome/horse_10c_norm.csv", ',', true);
CSVConfig const kOdTestNormIris = CreateCsvConfig("od_norm_data/metanome/iris_norm.csv", ',', true);
CSVConfig const kIndTestWide2 = CreateCsvConfig("ind_data/TestWide2.csv", ',', false);
CSVConfig const kIndTestEmpty = CreateCsvConfig("ind_data/Empty.csv", ',', true);
CSVConfig const kIndTestPlanets = CreateCsvConfig("ind_data/Planets.csv", ',', false);
CSVConfig const kIndTest3aryInds = CreateCsvConfig("ind_data/Test-3ary-inds.csv", ',', false);
CSVConfig const kIndTestTableFirst = CreateCsvConfig("ind_data/two_tables/first.csv", ',', false);
CSVConfig const kIndTestTableSecond = CreateCsvConfig("ind_data/two_tables/second.csv", ',', false);
CSVConfig const kIndTestNulls = CreateCsvConfig("INDTestNulls.csv", ',', true);
CSVConfig const kTestZeros = CreateCsvConfig("TestZeros.csv", ',', true);
CSVConfig const kNullEmpty = CreateCsvConfig("NullEmpty.csv", ',', true);
CSVConfig const kSimpleTypes = CreateCsvConfig("SimpleTypes.csv", ',', true);
CSVConfig const kRulesBook = CreateCsvConfig("transactional_data/rules-book.csv", ',', false);
CSVConfig const kRulesPresentationExtended =
        CreateCsvConfig("transactional_data/rules-presentation-extended.csv", ',', false);
CSVConfig const kRulesPresentation =
        CreateCsvConfig("transactional_data/rules-presentation.csv", ',', false);
CSVConfig const kRulesSynthetic2 =
        CreateCsvConfig("transactional_data/rules-synthetic-2.csv", ',', false);
CSVConfig const kRulesKaggleRows =
        CreateCsvConfig("transactional_data/rules-kaggle-rows.csv", ',', true);
CSVConfig const kTennis = CreateCsvConfig("cfd_data/tennis.csv", ',', true);
CSVConfig const kMushroom = CreateCsvConfig("cfd_data/mushroom.csv", ',', true);
CSVConfig const kTestDataStats = CreateCsvConfig("TestDataStats.csv", ',', false);
CSVConfig const kTestMetric = CreateCsvConfig("TestMetric.csv", ',', true);
CSVConfig const kBernoulliRelation = CreateCsvConfig("BernoulliRelation.csv", ',', true);
CSVConfig const kACShippingDates = CreateCsvConfig("ACShippingDates.csv", ',', true);
CSVConfig const kSimpleTypos = CreateCsvConfig("SimpleTypos.csv", ',', true);
CSVConfig const kTest1 = CreateCsvConfig("Test1.csv", ',', true);
CSVConfig const kProbeTest1 = CreateCsvConfig("ProbeTest1.csv", ',', true);
CSVConfig const kProbeTest2 = CreateCsvConfig("ProbeTest2.csv", ',', true);
CSVConfig const kTestParse = CreateCsvConfig("TestParse.csv", ',', false);
CSVConfig const kODnorm6 = CreateCsvConfig("OD_norm6.csv", ',', true);
CSVConfig const kTestDD = CreateCsvConfig("TestDD.csv", ',', true);
CSVConfig const kTestDD1 = CreateCsvConfig("TestDD1.csv", ',', true);
CSVConfig const kTestDD2 = CreateCsvConfig("TestDD2.csv", ',', true);
CSVConfig const kTestDD3 = CreateCsvConfig("TestDD3.csv", ',', true);
CSVConfig const kTestDif = CreateCsvConfig("dif_tables/TestDif.csv", ',', true);
CSVConfig const kTestDif1 = CreateCsvConfig("dif_tables/TestDif1.csv", ',', true);
CSVConfig const kTestDif2 = CreateCsvConfig("dif_tables/TestDif2.csv", ',', true);
CSVConfig const kTestDif3 = CreateCsvConfig("dif_tables/TestDif3.csv", ',', true);
CSVConfig const kSimpleTypes1 = CreateCsvConfig("SimpleTypes1.csv", ',', true);
}  // namespace tests
