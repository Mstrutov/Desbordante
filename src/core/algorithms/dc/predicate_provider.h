#pragma once

#include "dc/base_provider.h"
#include "dc/column_operand.h"
#include "operator.h"
#include "predicate.h"

namespace model {

/**
 * @brief Singleton storage for Predicate objects.
 *
 * A Predicate represents a relational condition between two rows from potentially different
 * columns, denoted as "t.A_i op s.A_j", where t and s represent different row identifiers, A_i
 * and A_j represent columns (which may be the same or different), and 'op' is a relational operator
 * (one of <, <=, >, >=, ==, !=). This class manages a centralized repository of Predicate objects.
 *
 * Usage:
 * To access or create a Predicate, use PredicateProvider::GetInstance()->GetPredicate(op, operand1,
 * operand2), where 'op' is the relational operator, and 'operand1' and 'operand2' are the column
 * operands involved.
 */
class PredicateProvider : public BaseProvider<PredicateProvider> {
private:
    using ColumnMap = std::unordered_map<ColumnOperand, Predicate>;
    using OperatorMap = std::unordered_map<ColumnOperand, ColumnMap>;
    // predicates_[op][col1][col2] corresponds to the related Predicate object
    std::unordered_map<Operator, OperatorMap> predicates_;

    friend BaseProvider<PredicateProvider>;

    static std::string ClassName() {
        return "PredicateProvider";
    }

    static void Clear() {
        instance_->predicates_.clear();
    }

public:
    PredicatePtr GetPredicate(Operator const& op, ColumnOperand const& left,
                              ColumnOperand const& right);
};

}  // namespace model
