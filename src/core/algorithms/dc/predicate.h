#pragma once

#include <span>

#include <boost/functional/hash.hpp>

#include "column_operand.h"
#include "operator.h"
#include "table/typed_column_data.h"

namespace model {

class Predicate;
using PredicatePtr = Predicate const*;
using PredicatesVector = std::vector<PredicatePtr>;
using PredicatesSpan = std::span<const PredicatePtr>;

/**
 * @brief Represents a predicate for FastADC.
 *
 * FastADC processes Denial Constraints (DCs) that involve comparisons between
 * pairs of rows within a dataset. A typical DC example, derived from a Functional
 * Dependency such as A -> B, is expressed as: `forall t, s in r, not (t.A = s.A and t.B != s.B)`
 * This denotes that for any pair of rows in the relation, it should not be the case
 * that while the values in column "A" are equal, the values in column "B" are unequal.
 *
 * A predicate in this context (e.g., t.A == s.A) comprises three elements to be fully
 * represented: the column operand from the first tuple ("t.A"), the comparison operator
 * ("="), and the column operand from the second tuple ("s.A").
 */
class Predicate {
private:
    Operator op_;
    ColumnOperand l_;
    ColumnOperand r_;

    mutable PredicatePtr symmetric_{};
    mutable PredicatePtr operator_symmetric_{};
    mutable PredicatePtr inv_TS_{};
    mutable PredicatePtr inverse_{};
    mutable std::vector<PredicatePtr> implications_;

public:
    Predicate(Operator const& op, ColumnOperand const& l, ColumnOperand const& r)
        : op_(op), l_(l), r_(r) {}

    // FIXME: mb pass some table representation other than vector of columns data?
    bool Satisfies(std::vector<model::TypedColumnData>& col_data, int t, int s) const;

    PredicatePtr GetSymmetric() const;

    PredicatePtr GetInvTS() const;

    PredicatePtr GetInverse() const;

    std::vector<PredicatePtr> GetImplications() const;

    Operator GetOperator() const {
        return op_;
    }

    ColumnOperand GetLeftOperand() const {
        return l_;
    }

    ColumnOperand GetRightOperand() const {
        return r_;
    }

    bool IsCrossColumn() const {
        return l_.GetColumn() != r_.GetColumn();
    }

    bool HasSameOperandsAs(Predicate const& rhs) const {
        return l_ == rhs.GetLeftOperand() && r_ == rhs.GetRightOperand();
    }

    bool HasSameOperandsAs(PredicatePtr rhs) const {
        return l_ == rhs->GetLeftOperand() && r_ == rhs->GetRightOperand();
    }

    bool operator==(Predicate const& rhs) const {
        return op_ == rhs.op_ && l_ == rhs.l_ && r_ == rhs.r_;
    }

    std::string ToString() const {
        return l_.ToString() + " " + op_.ToString() + " " + r_.ToString();
    }
};

/** Create predicate object and return pointer to it or obtain it from cache */
PredicatePtr GetPredicate(Operator const& op, ColumnOperand const& l, ColumnOperand const& r);

PredicatePtr GetPredicateByType(PredicatesSpan predicates, OperatorType type);

}  // namespace model

namespace std {
template <>
struct hash<model::Predicate> {
    size_t operator()(model::Predicate const& k) const noexcept {
        std::size_t seed = 0;
        boost::hash_combine(seed, k.GetOperator());
        boost::hash_combine(seed, k.GetLeftOperand());
        boost::hash_combine(seed, k.GetRightOperand());
        return seed;
    }
};

template <>
struct hash<std::shared_ptr<model::Predicate>> {
    size_t operator()(std::shared_ptr<model::Predicate> const& k) const {
        return std::hash<model::Predicate>()(*k);
    }
};
};  // namespace std
