#include "dc/predicate.h"

#include "dc/predicate_provider.h"

namespace model {

PredicatePtr GetPredicate(Operator const& op, ColumnOperand const& l, ColumnOperand const& r) {
    return PredicateProvider::GetInstance()->GetPredicate(op, l, r);
}

PredicatePtr GetPredicateByType(PredicatesSpan predicates, OperatorType type) {
    auto it = std::find_if(predicates.begin(), predicates.end(),
                           [type](PredicatePtr p) { return p->GetOperator() == type; });

    return it == predicates.end() ? nullptr : *it;
}

bool Predicate::Satisfies(std::vector<model::TypedColumnData>& col_data, int t, int s) const {
    TypedColumnData const& lhs = col_data[l_.GetColumn()->GetIndex()];
    TypedColumnData const& rhs = col_data[r_.GetColumn()->GetIndex()];

    // Assumes that types in both columns are the same (and they should)
    Type const& type = lhs.GetType();

    std::byte const* l_val = lhs.GetValue(l_.GetTuple() ? t : s);
    std::byte const* r_val = rhs.GetValue(r_.GetTuple() ? t : s);

    return op_.Eval(l_val, r_val, type);
}

PredicatePtr Predicate::GetSymmetric() const {
    if (!symmetric_) {
        symmetric_ = GetPredicate(op_.GetSymmetric(), r_, l_);
    }
    return symmetric_;
}

PredicatePtr Predicate::GetInvTS() const {
    if (!inv_TS_) {
        inv_TS_ = GetPredicate(op_, r_.GetInvTS(), l_.GetInvTS());
    }
    return inv_TS_;
}

PredicatePtr Predicate::GetInverse() const {
    if (!inverse_) {
        inverse_ = GetPredicate(op_.GetInverse(), l_, r_);
    }
    return inverse_;
}

std::vector<PredicatePtr> Predicate::GetImplications() const {
    if (implications_.empty()) {
        auto op_implications = op_.GetImplications();
        for (auto const& op_implication : op_implications) {
            implications_.push_back(GetPredicate(op_implication, l_, r_));
        }
    }
    return implications_;
}

}  // namespace model
