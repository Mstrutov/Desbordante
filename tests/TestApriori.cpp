#include <filesystem>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "EnumerationTree.h"
#include "AR.h"
#include "Datasets.h"

namespace fs = std::filesystem;

testing::AssertionResult checkFrequentListsEquality(std::list<std::set<std::string>> const& actual,
                                                    std::set<std::set<std::string>>  const& expected) {
    if (actual.size() != expected.size()) {
        return ::testing::AssertionFailure() << "count of frequent itemsets does not match: expected "
                                             << expected.size() << ", got " << actual.size();
    } else {
        for (auto const& itemset : actual) {
            if (expected.find(itemset) == expected.end()) {
                return testing::AssertionFailure() << "generated itemset not found in expected";
            }
        }
        return ::testing::AssertionSuccess();
    }
}

testing::AssertionResult checkAssociationRulesListsEquality(std::list<ARStrings> const& actual,
                                                            std::set<std::pair<std::set<std::string>,
                                                                               std::set<std::string>>> const& expected) {
    if (actual.size() != expected.size()) {
        return ::testing::AssertionFailure() << "count of generated rules does not match: expected "
                                             << expected.size() << ", got " << actual.size();
    } else {
        for (auto const& rule : actual) {
            std::set<std::string> lhs(rule.left.begin(), rule.left.end());
            std::set<std::string> rhs(rule.right.begin(), rule.right.end());

            if (expected.find(std::make_pair(std::move(lhs), std::move(rhs))) == expected.end()) {
                return testing::AssertionFailure() << "generated rule does not present in expected";
            }
        }
        return ::testing::AssertionSuccess();
    }
}

class ARAlgorithmTest : public ::testing::Test {
protected:
    std::unique_ptr<ARAlgorithm> createAlgorithmInstance(
            double minsup, double minconf,
            std::filesystem::path const& path,
            std::shared_ptr<InputFormat> inputFormat,
            char separator = ',', bool hasHeader = true) {
        ARAlgorithm::Config const config = {path, separator, hasHeader, std::move(inputFormat), minsup, minconf};
        return std::make_unique<EnumerationTree>(config);
    }
};

TEST_F(ARAlgorithmTest, BookDataset) {
    auto const path = fs::current_path() / "inputData" / "transactionalData" / "rules-book.csv";
    auto inputParams = std::make_shared<Singular>(0, 1);
    auto algorithm = createAlgorithmInstance(0.3, 0.5, path, std::move(inputParams), ',', false);
    algorithm->Execute();
    auto const actualFrequent = algorithm->getAllFrequent();
    std::set<std::set<std::string>> const expectedFrequent = {
             {"Bread"}, {"Milk"}, {"Eggs"}, {"Cheese"}, {"Yogurt"},
             {"Bread", "Milk"}, {"Eggs", "Milk"}, {"Cheese", "Milk"}, {"Eggs", "Yogurt"}, {"Milk", "Yogurt"},
             {"Eggs", "Milk", "Yogurt"}
    };
    ASSERT_TRUE(checkFrequentListsEquality(actualFrequent, expectedFrequent));

    auto const actualRules = algorithm->arList();
    std::set<std::pair<std::set<std::string>, std::set<std::string>>> expectedRules = {
            {{"Bread"}, {"Milk"}},
            {{"Cheese"}, {"Milk"}},
            {{"Milk"}, {"Eggs"}},
            {{"Eggs"}, {"Milk"}},
            {{"Yogurt"}, {"Eggs"}},
            {{"Eggs"}, {"Yogurt"}},
            {{"Yogurt"}, {"Milk"}},
            {{"Milk"}, {"Yogurt"}},
            {{"Yogurt", "Milk"}, {"Eggs"}},
            {{"Yogurt", "Eggs"}, {"Milk"}},
            {{"Milk", "Eggs"}, {"Yogurt"}},
            {{"Yogurt"}, {"Milk", "Eggs"}},
            {{"Eggs"}, {"Yogurt", "Milk"}}
    };
    ASSERT_TRUE(checkAssociationRulesListsEquality(actualRules, expectedRules));
}

TEST_F(ARAlgorithmTest, PresentationExtendedDataset) {
    auto const path = fs::current_path() / "inputData" / "transactionalData" / "rules-presentation-extended.csv";
    auto inputParams = std::make_shared<Singular>(0, 1);
    auto algorithm = createAlgorithmInstance(0.6, 0, path, std::move(inputParams), ',', false);
    algorithm->Execute();
    auto const actual = algorithm->getAllFrequent();
    std::set<std::set<std::string>> const expected = {
            {"Bread"}, {"Milk"}, {"Diaper"}, {"Beer"},
            {"Bread", "Milk"}, {"Diaper", "Beer"}, {"Beer", "Milk"}, {"Bread", "Beer"}, {"Milk", "Diaper"}, {"Bread", "Diaper"},
            {"Bread", "Diaper", "Beer"}, {"Milk", "Diaper", "Beer"}
    };
    ASSERT_TRUE(checkFrequentListsEquality(actual, expected));
}

TEST_F(ARAlgorithmTest, PresentationDataset) {
    auto const path = fs::current_path() / "inputData" / "transactionalData" / "rules-presentation.csv";
    auto inputParams = std::make_shared<Singular>(0, 1);
    auto algorithm = createAlgorithmInstance(0.6, 0, path, std::move(inputParams), ',', false);
    algorithm->Execute();

    auto const actual = algorithm->getAllFrequent();
    std::set<std::set<std::string>> const expected = {
            {"Bread"}, {"Milk"}, {"Diaper"}, {"Beer"},
            {"Bread", "Milk"}, {"Diaper", "Beer"}, {"Milk", "Diaper"}, {"Bread", "Diaper"}
    };
    ASSERT_TRUE(checkFrequentListsEquality(actual, expected));

    auto const actualRules = algorithm->arList();
    std::set<std::pair<std::set<std::string>, std::set<std::string>>> expectedRules = {
            {{"Bread"}, {"Milk"}},
            {{"Milk"}, {"Bread"}},
            {{"Diaper"}, {"Beer"}},
            {{"Beer"}, {"Diaper"}},
            {{"Diaper"}, {"Milk"}},
            {{"Milk"}, {"Diaper"}},
            {{"Diaper"}, {"Bread"}},
            {{"Bread"}, {"Diaper"}}
    };
    ASSERT_TRUE(checkAssociationRulesListsEquality(actualRules, expectedRules))
            << "conf=0: generated not all of the combinations from the frequent itemsets";
}

TEST_F(ARAlgorithmTest, SynteticDatasetWithPruning) {
    auto const path = fs::current_path() / "inputData" / "transactionalData" / "rules-synthetic-2.csv";
    auto inputParams = std::make_shared<Singular>(0, 1);
    auto algorithm = createAlgorithmInstance(0.13, 1.00001, path, std::move(inputParams), ',', false);
    algorithm->Execute();

    auto const actual = algorithm->getAllFrequent();
    std::set<std::set<std::string>> const expected = {
            {"a"}, {"b"}, {"c"}, {"d"}, {"e"}, {"f"},
            {"a", "b"}, {"a", "c"}, {"a", "d"}, {"a", "f"}, {"b", "c"}, {"c", "d"}, {"c", "f"}, {"d", "f"},
            {"a", "c", "d"}, {"a", "c", "f"}, {"a", "d", "f"}, {"c", "d", "f"},
            {"a", "c", "d", "f"}
    };
    ASSERT_TRUE(checkFrequentListsEquality(actual, expected));

    auto const actualRules = algorithm->arList();
    std::set<std::pair<std::set<std::string>, std::set<std::string>>> expectedRules = {};
    ASSERT_TRUE(checkAssociationRulesListsEquality(actualRules, expectedRules))
             << "conf=1: generated some rules with the confidence value above one";
}

TEST_F(ARAlgorithmTest, KaggleDatasetWithTIDandHeader) {
    auto const path = fs::current_path() / "inputData" / "transactionalData" / "rules-kaggle-rows.csv";
    auto inputParams = std::make_shared<Tabular>(true);
    auto algorithm = createAlgorithmInstance(0.1, 0.5, path, std::move(inputParams),',', true);
    algorithm->Execute();

    auto const actualFrequent = algorithm->getAllFrequent();
    std::set<std::set<std::string>> const expectedFrequent = {
            {"BISCUIT"}, {"BOURNVITA"}, {"BREAD"}, {"COCK"}, {"COFFEE"}, {"CORNFLAKES"}, {"JAM"}, {"MAGGI"},
                {"MILK"}, {"SUGER"}, {"TEA"},
            {"BISCUIT", "BREAD"}, {"COCK", "BISCUIT"}, {"BISCUIT", "COFFEE"}, {"CORNFLAKES", "BISCUIT"}, {"MAGGI", "BISCUIT"},
                {"MILK", "BISCUIT"}, {"BISCUIT", "TEA"}, {"BREAD", "BOURNVITA"}, {"SUGER", "BOURNVITA"}, {"TEA", "BOURNVITA"},
                {"COFFEE", "BREAD"}, {"JAM", "BREAD"}, {"MAGGI", "BREAD"}, {"MILK", "BREAD"}, {"BREAD", "SUGER"}, {"TEA", "BREAD"},
                {"COCK", "COFFEE"}, {"CORNFLAKES", "COCK"}, {"CORNFLAKES", "COFFEE"}, {"COFFEE", "SUGER"}, {"CORNFLAKES", "MILK"},
                {"CORNFLAKES", "TEA"}, {"JAM", "MAGGI"}, {"MAGGI", "TEA"},
            {"MILK", "BISCUIT", "BREAD"}, {"COCK", "BISCUIT", "COFFEE"}, {"CORNFLAKES", "COCK", "BISCUIT"},
                {"CORNFLAKES", "BISCUIT", "COFFEE"}, {"MAGGI", "BISCUIT", "TEA"}, {"TEA", "BREAD", "BOURNVITA"},
                {"BREAD", "COFFEE", "SUGER"}, {"JAM", "MAGGI", "BREAD"}, {"MAGGI", "TEA", "BREAD"},
                {"CORNFLAKES", "COCK", "COFFEE"},
            {"CORNFLAKES", "COCK", "BISCUIT", "COFFEE"},
    };
    ASSERT_TRUE(checkFrequentListsEquality(actualFrequent, expectedFrequent));

    auto const actualRules = algorithm->arList();
    std::set<std::pair<std::set<std::string>, std::set<std::string>>> expectedRules = {
            {{"BISCUIT"}, {"BREAD"}},
            {{"COCK"}, {"BISCUIT"}},
            {{"CORNFLAKES"}, {"BISCUIT"}},
            {{"BOURNVITA"}, {"BREAD"}},
            {{"BOURNVITA"}, {"SUGER"}},
            {{"BOURNVITA"}, {"TEA"}},
            {{"JAM"}, {"BREAD"}},
            {{"MAGGI"}, {"BREAD"}},
            {{"MILK"}, {"BREAD"}},
            {{"SUGER"}, {"BREAD"}},
            {{"TEA"}, {"BREAD"}},
            {{"COCK"}, {"COFFEE"}},
            {{"COCK"}, {"CORNFLAKES"}},
            {{"COFFEE"}, {"CORNFLAKES"}},
            {{"CORNFLAKES"}, {"COFFEE"}},
            {{"COFFEE"}, {"SUGER"}},
            {{"SUGER"}, {"COFFEE"}},
            {{"JAM"}, {"MAGGI"}},
            {{"MAGGI"}, {"TEA"}},
            {{"TEA"}, {"MAGGI"}},
            {{"MILK", "BISCUIT"}, {"BREAD"}},
            {{"MILK", "BREAD"}, {"BISCUIT"}},
            {{"BISCUIT", "BREAD"}, {"MILK"}},
            {{"COCK", "BISCUIT"}, {"COFFEE"}},
            {{"COCK", "COFFEE"}, {"BISCUIT"}},
            {{"BISCUIT", "COFFEE"}, {"COCK"}},
            {{"COCK"}, {"BISCUIT", "COFFEE"}},
            {{"COCK", "BISCUIT"}, {"CORNFLAKES"}},
            {{"COCK", "CORNFLAKES"}, {"BISCUIT"}},
            {{"BISCUIT", "CORNFLAKES"}, {"COCK"}},
            {{"COCK"}, {"BISCUIT", "CORNFLAKES"}},
            {{"BISCUIT", "COFFEE"}, {"CORNFLAKES"}},
            {{"BISCUIT", "CORNFLAKES"}, {"COFFEE"}},
            {{"COFFEE", "CORNFLAKES"}, {"BISCUIT"}},
            {{"BISCUIT", "MAGGI"}, {"TEA"}},
            {{"BISCUIT", "TEA"}, {"MAGGI"}},
            {{"MAGGI", "TEA"}, {"BISCUIT"}},
            {{"TEA", "BREAD"}, {"BOURNVITA"}},
            {{"BOURNVITA", "BREAD"}, {"TEA"}},
            {{"TEA", "BOURNVITA"}, {"BREAD"}},
            {{"BOURNVITA"}, {"TEA", "BREAD"}},
            {{"COFFEE", "SUGER"}, {"BREAD"}},
            {{"COFFEE", "BREAD"}, {"SUGER"}},
            {{"SUGER", "BREAD"}, {"COFFEE"}},
            {{"JAM", "MAGGI"}, {"BREAD"}},
            {{"JAM", "BREAD"}, {"MAGGI"}},
            {{"MAGGI", "BREAD"}, {"JAM"}},
            {{"JAM"}, {"MAGGI", "BREAD"}},
            {{"MAGGI", "TEA"}, {"BREAD"}},
            {{"MAGGI", "BREAD"}, {"TEA"}},
            {{"TEA", "BREAD"}, {"MAGGI"}},
            {{"COCK", "COFFEE"}, {"CORNFLAKES"}},
            {{"COCK", "CORNFLAKES"}, {"COFFEE"}},
            {{"COFFEE", "CORNFLAKES"}, {"COCK"}},
            {{"COCK"}, {"COFFEE", "CORNFLAKES"}},
            {{"COCK", "BISCUIT", "COFFEE"}, {"CORNFLAKES"}},
            {{"COCK", "BISCUIT", "CORNFLAKES"}, {"COFFEE"}},
            {{"COCK", "COFFEE", "CORNFLAKES"}, {"BISCUIT"}},
            {{"BISCUIT", "COFFEE", "CORNFLAKES"}, {"COCK"}},
            {{"COCK", "BISCUIT"}, {"COFFEE", "CORNFLAKES"}},
            {{"COCK", "COFFEE"}, {"BISCUIT", "CORNFLAKES"}},
            {{"COCK", "CORNFLAKES"}, {"BISCUIT", "COFFEE"}},
            {{"BISCUIT", "COFFEE"}, {"COCK", "CORNFLAKES"}},
            {{"BISCUIT", "CORNFLAKES"}, {"COCK", "COFFEE"}},
            {{"COFFEE", "CORNFLAKES"}, {"COCK", "BISCUIT"}},
            {{"COCK"}, {"BISCUIT", "COFFEE", "CORNFLAKES"}}
    };
    ASSERT_TRUE(checkAssociationRulesListsEquality(actualRules, expectedRules));
}