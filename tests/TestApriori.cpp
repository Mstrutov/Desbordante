#include <filesystem>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "EnumerationTree.h"
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

class ARAlgorithmTest : public ::testing::Test {
protected:
    std::unique_ptr<ARAlgorithm> createAlgorithmInstance(
            double minsup, double minconf,
            std::filesystem::path const& path,
            TransactionalInputFormat inputFormat = TransactionalInputFormat::TwoColumns,
            bool hasTID = false, char separator = ',', bool hasHeader = true) {
        return std::make_unique<EnumerationTree>(minsup, minconf, path, inputFormat, hasTID, separator, hasHeader);
    }
};

TEST_F(ARAlgorithmTest, BookDataset) {
    auto const path = fs::current_path() / "inputData" / "transactionalData" / "rules-book.csv";
    auto algorithm = createAlgorithmInstance(0.3, 0.4, path, TransactionalInputFormat::TwoColumns, false, ',', false);
    algorithm->execute();
    auto const actual = algorithm->getAllFrequent();
    std::set<std::set<std::string>> const expected = {
             {"Bread"}, {"Milk"}, {"Eggs"}, {"Cheese"}, {"Yogurt"},
             {"Bread", "Milk"}, {"Eggs", "Milk"}, {"Cheese", "Milk"}, {"Eggs", "Yogurt"}, {"Milk", "Yogurt"},
             {"Eggs", "Milk", "Yogurt"}
    };
    ASSERT_TRUE(checkFrequentListsEquality(actual, expected));
}

TEST_F(ARAlgorithmTest, PresentationExtendedDataset) {
    auto const path = fs::current_path() / "inputData" / "transactionalData" / "rules-presentation-extended.csv";
    auto algorithm = createAlgorithmInstance(0.6, 0, path, TransactionalInputFormat::TwoColumns, false, ',', false);
    algorithm->execute();
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
    auto algorithm = createAlgorithmInstance(0.6, 0, path, TransactionalInputFormat::TwoColumns, false, ',', false);
    algorithm->execute();
    auto const actual = algorithm->getAllFrequent();
    std::set<std::set<std::string>> const expected = {
            {"Bread"}, {"Milk"}, {"Diaper"}, {"Beer"},
            {"Bread", "Milk"}, {"Diaper", "Beer"}, {"Milk", "Diaper"}, {"Bread", "Diaper"}
    };
    ASSERT_TRUE(checkFrequentListsEquality(actual, expected));
}

TEST_F(ARAlgorithmTest, SynteticDatasetWithPruning) {
    auto const path = fs::current_path() / "inputData" / "transactionalData" / "rules-synthetic-2.csv";
    auto algorithm = createAlgorithmInstance(0.13, 0, path, TransactionalInputFormat::TwoColumns, false, ',', false);
    algorithm->execute();
    auto const actual = algorithm->getAllFrequent();
    std::set<std::set<std::string>> const expected = {
            {"a"}, {"b"}, {"c"}, {"d"}, {"e"}, {"f"},
            {"a", "b"}, {"a", "c"}, {"a", "d"}, {"a", "f"}, {"b", "c"}, {"c", "d"}, {"c", "f"}, {"d", "f"},
            {"a", "c", "d"}, {"a", "c", "f"}, {"a", "d", "f"}, {"c", "d", "f"},
            {"a", "c", "d", "f"}
    };
    ASSERT_TRUE(checkFrequentListsEquality(actual, expected));
}

TEST_F(ARAlgorithmTest, KaggleDatasetWithTIDandHeader) {
    auto const path = fs::current_path() / "inputData" / "transactionalData" / "rules-kaggle-rows.csv";
    auto algorithm = createAlgorithmInstance(0.1, 0.5, path, TransactionalInputFormat::ItemsetRows, true, ',', true);
    algorithm->execute();
    auto const actual = algorithm->getAllFrequent();
    std::set<std::set<std::string>> const expected = {
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
    ASSERT_TRUE(checkFrequentListsEquality(actual, expected));
}