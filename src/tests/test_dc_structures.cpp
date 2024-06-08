#include <memory>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "all_csv_configs.h"
#include "csv_parser/csv_parser.h"
#include "dc/column_operand.h"
#include "dc/operator.h"
#include "dc/predicate.h"
#include "dc/predicate_builder.h"
#include "table/column_layout_typed_relation_data.h"

namespace tests {

namespace fs = std::filesystem;
namespace mo = model;

class TestOperatorInt : public ::testing::Test {
private:
    std::unique_ptr<std::byte[]> left_owner_;
    std::unique_ptr<std::byte[]> right_owner_;

protected:
    std::unique_ptr<mo::IntType> type_;
    std::byte* left_ptr_;
    std::byte* right_ptr_;
    std::array<mo::Operator, 6> all_operators_ = {
            {mo::Operator(mo::OperatorType::kEqual), mo::Operator(mo::OperatorType::kUnequal),
             mo::Operator(mo::OperatorType::kGreater), mo::Operator(mo::OperatorType::kLess),
             mo::Operator(mo::OperatorType::kGreaterEqual),
             mo::Operator(mo::OperatorType::kLessEqual)}};

    void SetUp() override {
        type_ = std::make_unique<mo::IntType>();
        left_owner_.reset(type_->MakeValue(0));
        right_owner_.reset(type_->MakeValue(0));
        left_ptr_ = left_owner_.get();
        right_ptr_ = right_owner_.get();
    }

    void SetVals(int left, int right) {
        mo::Type::GetValue<int>(left_ptr_) = left;
        mo::Type::GetValue<int>(right_ptr_) = right;
    }
};

TEST_F(TestOperatorInt, Equal) {
    mo::Operator op(mo::OperatorType::kEqual);

    auto test = [&](int l, int r, bool res) {
        this->SetVals(l, r);

        EXPECT_EQ(res, op.Eval(this->left_ptr_, this->right_ptr_, *(this->type_)));
    };

    test(1, 1, true);
    test(1, 2, false);
}

TEST_F(TestOperatorInt, Greater) {
    mo::Operator op(mo::OperatorType::kGreater);

    auto test = [&](int l, int r, bool res) {
        this->SetVals(l, r);

        EXPECT_EQ(res, op.Eval(this->left_ptr_, this->right_ptr_, *(this->type_)));
    };

    test(1, 1, false);
    test(1, 2, false);
    test(-10, -20, true);
}

TEST_F(TestOperatorInt, Implication) {
    auto test = [&](int l, int r) {
        this->SetVals(l, r);

        for (auto& op : all_operators_) {
            auto implications = op.GetImplications();
            bool op_result = op.Eval(this->left_ptr_, this->right_ptr_, *(this->type_));

            for (auto& imp_op : implications) {
                if (op_result) {
                    EXPECT_TRUE(imp_op.Eval(this->left_ptr_, this->right_ptr_, *(this->type_)));
                }
            }
        }
    };

    test(0, 0);
    test(5, 4);
    test(5, -5);
    test(5, 5);
    test(10, 4);
    test(12313, -10);
}

TEST_F(TestOperatorInt, Inverse) {
    auto test = [&](int l, int r) {
        this->SetVals(l, r);

        for (auto& op : all_operators_) {
            auto inverse_op = op.GetInverse();
            bool op_result = op.Eval(this->left_ptr_, this->right_ptr_, *(this->type_));
            bool inverse_result =
                    inverse_op.Eval(this->left_ptr_, this->right_ptr_, *(this->type_));

            // The result of an operation and its inverse should be opposite
            EXPECT_NE(op_result, inverse_result);
        }
    };

    test(0, 0);
    test(5, 4);
    test(5, -5);
    test(5, 5);
    test(10, 4);
    test(12313, -10);
}

TEST_F(TestOperatorInt, Symmetric) {
    auto test = [&](int l, int r) {
        this->SetVals(l, r);

        for (auto& op : all_operators_) {
            auto symmetric_op = op.GetSymmetric();
            bool op_result = op.Eval(this->left_ptr_, this->right_ptr_, *(this->type_));

            bool symmetric_result =
                    symmetric_op.Eval(this->right_ptr_, this->left_ptr_, *(this->type_));

            EXPECT_EQ(op_result, symmetric_result);
        }
    };

    test(0, 0);
    test(5, 4);
    test(5, -5);
    test(5, 5);
    test(10, 4);
    test(12313, -10);
}

TEST_F(TestOperatorInt, Transitives) {
    auto test = [&](int a, int b, int c) {
        for (auto& op : all_operators_) {
            this->SetVals(a, b);
            bool ab_result = op.Eval(this->left_ptr_, this->right_ptr_, *(this->type_));

            auto transitives = op.GetTransitives();
            for (auto& trans_op : transitives) {
                this->SetVals(b, c);
                bool bc_trans_result =
                        trans_op.Eval(this->left_ptr_, this->right_ptr_, *(this->type_));

                if (ab_result && bc_trans_result) {
                    this->SetVals(a, c);
                    // If 'a op b' and 'b trans_op c', then 'a op c' should also be true
                    EXPECT_TRUE(op.Eval(this->left_ptr_, this->right_ptr_, *(this->type_)));
                }
            }
        }
    };

    test(1, 2, 3);
    test(3, 2, 1);
    test(2, 2, 2);
    test(-2, 0, 2);
    test(-3, -2, 1);
    test(-1, -5, -4);
    test(2, 2, 3);
    test(-1, -1, 0);
    test(1, 3, 3);
    test(-4, -2, -2);
    test(5, 3, 4);
    test(0, -1, 1);
    test(4, 2, 2);
    test(1, 0, 0);
    test(1000, 2000, 3000);
    test(-1000, -500, 0);
    test(10000, -5000, 5000);
    test(3, 3, 2);
    test(-1, -1, -2);
    test(2, 3, 1);
    test(-3, -1, -4);
    test(3, 1, 2);
    test(0, -2, -1);
}

TEST(TestOperatorString, Compare) {
    std::unique_ptr<mo::Type> type(mo::CreateType(mo::TypeId::kString, true));
    mo::StringType const* s = static_cast<mo::StringType const*>(type.get());

    mo::Operator eq(mo::OperatorType::kEqual);
    mo::Operator neq(mo::OperatorType::kUnequal);

    auto test = [&](std::string l, std::string r, bool res, mo::Operator& op) {
        using Owner = std::unique_ptr<std::byte[], mo::StringTypeDeleter>;

        Owner l_owner(s->MakeValue(l), s->GetDeleter());
        Owner r_owner(s->MakeValue(r), s->GetDeleter());

        std::byte* l_val = l_owner.get();
        std::byte* r_val = r_owner.get();

        EXPECT_EQ(res, op.Eval(l_val, r_val, *s));
    };

    test("sfsdf", "sfsdf", true, eq);
    test("sfsdf", "sfsdf", false, neq);
    test("abc", "cba", false, eq);
    test("abc", "cba", true, neq);
}

TEST(Predicate, PredicateCreatesCorrectly) {
    CSVParser parser{kTmpDC};
    std::unique_ptr<model::ColumnLayoutTypedRelationData> table =
            mo::ColumnLayoutTypedRelationData::CreateFrom(parser, true);
    std::vector<mo::TypedColumnData> col_data = std::move(table->GetColumnData());
    Column const *first = col_data[0].GetColumn(), *second = col_data[1].GetColumn();
    // FIXME: temprorary to make test work
    model::PredicateBuilder pbuilder;

    mo::PredicatePtr s_a_less_t_b =
            mo::GetPredicate(mo::Operator(mo::OperatorType::kLess), mo::ColumnOperand(first, true),
                             mo::ColumnOperand(second, false));

    EXPECT_TRUE(s_a_less_t_b->Satisfies(col_data, 0, 1));
    EXPECT_TRUE(s_a_less_t_b->Satisfies(col_data, 1, 0));

    mo::PredicatePtr s_a_neq_t_a =
            GetPredicate(mo::Operator(mo::OperatorType::kUnequal), mo::ColumnOperand(first, true),
                         mo::ColumnOperand(first, false));

    EXPECT_FALSE(s_a_neq_t_a->Satisfies(col_data, 0, 1));
    EXPECT_FALSE(s_a_neq_t_a->Satisfies(col_data, 1, 0));
}

}  // namespace tests
