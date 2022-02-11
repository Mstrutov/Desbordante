#include <iostream>
#include <thread>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ColumnLayoutRelationData.h"
#include "ListAgreeSetSample.h"
#include "IdentifierSet.h"
#include "AgreeSetFactory.h"

using std::deque, std::vector, std::cout, std::endl, std::unique_ptr, util::AgreeSetFactory,
    util::MCGenMethod, util::AgreeSetsGenMethod;
using ::testing::ContainerEq, ::testing::Eq;

namespace fs = std::filesystem;

TEST(pliChecker, first) {
    deque<vector<int>> ans = {
        {0, 2, 8, 11}, {1, 5, 9}, {4, 14}, {6, 7, 18}, {10, 17}  // null
    };

    auto path = fs::current_path().append("inputData").append("Test1.csv");
    deque<vector<int>> index;
    try {
        CSVParser csv_parser(path);
        auto test = ColumnLayoutRelationData::CreateFrom(csv_parser, true);
        auto column_data = test->GetColumnData(0);
        index = column_data.GetPositionListIndex()->GetIndex();
    }
    catch (std::runtime_error& e) {
        cout << "Exception raised in test: " << e.what() << endl;
        FAIL();
    }
    ASSERT_THAT(index, ContainerEq(ans));
}

TEST(pliChecker, second) {
    deque<vector<int>> ans = {
        {0, 2, 8, 11},
        {1, 5, 9},
        {4, 14},
        {6, 7, 18},
    };
    deque<vector<int>> index;
    try {
        auto path = fs::current_path().append("inputData").append("Test1.csv");
        CSVParser csv_parser(path);
        auto test = ColumnLayoutRelationData::CreateFrom(csv_parser, false);
        auto column_data = test->GetColumnData(0);
        index = column_data.GetPositionListIndex()->GetIndex();
    }
    catch (std::runtime_error& e) {
        cout << "Exception raised in test: " << e.what() << endl;
        FAIL();
    }
    ASSERT_THAT(index, ContainerEq(ans));
}

TEST(pliIntersectChecker, first) {
    deque<vector<int>> ans = {{2, 5}};
    std::shared_ptr<util::PositionListIndex> intersection;

    try {
        auto path = fs::current_path().append("inputData");
        CSVParser csv_parser_1(path / "ProbeTest1.csv");
        CSVParser csv_parser_2(path / "ProbeTest2.csv");

        auto test1 = ColumnLayoutRelationData::CreateFrom(csv_parser_1, false);
        auto test2 = ColumnLayoutRelationData::CreateFrom(csv_parser_2, false);
        auto pli1 = test1->GetColumnData(0).GetPositionListIndex();
        auto pli2 = test2->GetColumnData(0).GetPositionListIndex();

        intersection = pli1->Intersect(pli2);
    }
    catch (std::runtime_error& e) {
        cout << "Exception raised in test: " << e.what() << endl;
        FAIL();
    }
    ASSERT_THAT(intersection->GetIndex(), ContainerEq(ans));
}

TEST(testingBitsetToLonglong, first) {
    size_t encoded_num = 1254;
    boost::dynamic_bitset<> simple_bitset{20, encoded_num};

    auto res_vector = *util::ListAgreeSetSample::BitSetToLongLongVector(simple_bitset);
    ASSERT_EQ(res_vector.size(), 1);
    for (auto long_long_repr : res_vector)
        ASSERT_EQ(encoded_num, long_long_repr);
}

TEST(IdentifierSetTest, Computation) {
    std::set<std::string> id_sets;
    std::set<std::string> id_sets_ans = {
        "[(A, 0), (B, 1), (C, 1), (D, 1), (E, 1), (F, 1)]",
        "[(A, 1), (B, 1), (C, 1), (D, 2), (E, 0), (F, 1)]",
        "[(A, 1), (B, 1), (C, 2), (D, 1), (E, 1), (F, 0)]",
        "[(A, 1), (B, 2), (C, 2), (D, 1), (E, 1), (F, 1)]",
        "[(A, 1), (B, 2), (C, 1), (D, 2), (E, 1), (F, 1)]",
        "[(A, 1), (B, 2), (C, 2), (D, 2), (E, 1), (F, 1)]"
    };

    try {
        auto path = fs::current_path().append("inputData").append("BernoulliRelation.csv");
        CSVParser parser(path);
        auto relation = ColumnLayoutRelationData::CreateFrom(parser, false);

        for (unsigned i = 0; i < relation->GetNumRows(); ++i) {
            id_sets.insert(util::IdentifierSet(relation.get(), i).ToString());
        }
    }
    catch (std::runtime_error const& e) {
        cout << "Exception raised in test: " << e.what() << endl;
        FAIL();
    }
    ASSERT_THAT(id_sets_ans, ContainerEq(id_sets));
}

TEST(IdentifierSetTest, Intersection) {
    std::set<std::string> intersection_actual; // id set intersection result
    std::set<std::string> intersection_ans = {
        "[A D F]",
        "[A B]",
        "[D E F]",
        "[A E]",
        "[C E F]",
        "[E F]",
        "[B C F]",
        "[A B E F]",
        "[A C D E]",
        "[B D E]",
        "[A B C E F]",
        "[A F]",
        "[A B D E F]",
        "[A C D F]",
        "[A C E]"
    };

    try {
        auto path = fs::current_path().append("inputData").append("BernoulliRelation.csv");
        CSVParser parser(path);
        auto relation = ColumnLayoutRelationData::CreateFrom(parser, false);
        std::vector<util::IdentifierSet> id_sets;

        for (unsigned i = 0; i < relation->GetNumRows(); ++i) {
            id_sets.emplace_back(relation.get(), i);
        }

        auto back_it = std::prev(id_sets.end());
        for (auto p = id_sets.begin(); p != back_it; ++p) {
            for (auto q = std::next(p); q != id_sets.end(); ++q) {
                intersection_actual.insert(p->Intersect(*q).ToString());
            }
        }
    }
    catch (std::runtime_error const& e) {
        cout << "Exception raised in test: " << e.what() << endl;
        FAIL();
    }
    ASSERT_THAT(intersection_ans, ContainerEq(intersection_actual));
}

void TestAgreeSetFactory(AgreeSetFactory::Configuration c) {
    std::set<std::string> agree_sets_actual; // id set intersection result
    std::set<std::string> agree_sets_ans = {
        "[A D F]",
        "[A B]",
        "[D E F]",
        "[A E]",
        "[C E F]",
        "[E F]",
        "[B C F]",
        "[A B E F]",
        "[]",
        "[A C D E]",
        "[B D E]",
        "[A B C E F]",
        "[A F]",
        "[A B D E F]",
        "[A C D F]",
        "[A C E]"
    };

    try {
        auto path = fs::current_path().append("inputData").append("BernoulliRelation.csv");
        CSVParser parser(path);
        auto relation = ColumnLayoutRelationData::CreateFrom(parser, false);
        AgreeSetFactory factory(relation.get(), c);
        for (util::AgreeSet const& agree_set : factory.GenAgreeSets()) {
            agree_sets_actual.insert(agree_set.ToString());
        }
    }
    catch (std::runtime_error const& e) {
        cout << "Exception raised in test: " << e.what() << endl;
        FAIL();
    }
    ASSERT_THAT(agree_sets_ans, ContainerEq(agree_sets_actual));
}

TEST(AgreeSetFactoryTest, UsingVectorOfIDSets) {
    AgreeSetFactory::Configuration c(AgreeSetsGenMethod::kUsingVectorOfIDSets);
    TestAgreeSetFactory(c);
}

TEST(AgreeSetFactoryTest, UsingMapOfIDSets) {
    AgreeSetFactory::Configuration c(AgreeSetsGenMethod::kUsingMapOfIDSets);
    TestAgreeSetFactory(c);
}

TEST(AgreeSetFactoryTest, UsingGetAgreeSet) {
    AgreeSetFactory::Configuration c(AgreeSetsGenMethod::kUsingGetAgreeSet);
    TestAgreeSetFactory(c);
}

TEST(AgreeSetFactoryTest, UsingMCAndGetAgreeSet) {
    AgreeSetFactory::Configuration c(AgreeSetsGenMethod::kUsingMCAndGetAgreeSet);
    TestAgreeSetFactory(c);
}

TEST(AgreeSetFactoryTest, UsingHandleEqvClass) {
    AgreeSetFactory::Configuration c(MCGenMethod::kUsingHandleEqvClass);
    TestAgreeSetFactory(c);
}

TEST(AgreeSetFactoryTest, UsingCalculateSupersets) {
    AgreeSetFactory::Configuration c(MCGenMethod::kUsingCalculateSupersets);
    TestAgreeSetFactory(c);
}

TEST(AgreeSetFactoryTest, UsingHandlePartition) {
    AgreeSetFactory::Configuration c(MCGenMethod::kUsingHandlePartition);
    TestAgreeSetFactory(c);
}

#if 0
TEST(AgreeSetFactoryTest, MCGenParallel) {
    AgreeSetFactory::Configuration c(AgreeSetsGenMethod::kUsingVectorOfIDSets,
                                     MCGenMethod::kParallel,
                                     std::thread::hardware_concurrency());
    testAgreeSetFactory(c);
}
#endif

