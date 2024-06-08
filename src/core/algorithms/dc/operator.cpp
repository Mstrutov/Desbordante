#include "algorithms/dc/operator.h"

namespace model {

bool Operator::Eval(std::byte const* v1, std::byte const* v2, Type const& type) const {
    CompareResult cr = type.Compare(v1, v2);

    switch (op_) {
        case OperatorType::kEqual:
            return cr == CompareResult::kEqual;
        case OperatorType::kUnequal:
            return cr != CompareResult::kEqual;
        case OperatorType::kGreater:
            return cr == CompareResult::kGreater;
        case OperatorType::kLess:
            return cr == CompareResult::kLess;
        case OperatorType::kGreaterEqual:
            return cr == CompareResult::kGreater || cr == CompareResult::kEqual;
        case OperatorType::kLessEqual:
            return cr == CompareResult::kLess || cr == CompareResult::kEqual;
        default:
            return false;
    }
}

std::unordered_map<OperatorType, OperatorType> const Operator::InitializeInverseMap() {
    return {{OperatorType::kEqual, OperatorType::kUnequal},
            {OperatorType::kUnequal, OperatorType::kEqual},
            {OperatorType::kGreater, OperatorType::kLessEqual},
            {OperatorType::kLess, OperatorType::kGreaterEqual},
            {OperatorType::kGreaterEqual, OperatorType::kLess},
            {OperatorType::kLessEqual, OperatorType::kGreater}};
}

std::unordered_map<OperatorType, OperatorType> const Operator::InitializeSymmetricMap() {
    return {{OperatorType::kEqual, OperatorType::kEqual},
            {OperatorType::kUnequal, OperatorType::kUnequal},
            {OperatorType::kGreater, OperatorType::kLess},
            {OperatorType::kLess, OperatorType::kGreater},
            {OperatorType::kGreaterEqual, OperatorType::kLessEqual},
            {OperatorType::kLessEqual, OperatorType::kGreaterEqual}};
}

std::unordered_map<OperatorType, std::vector<OperatorType>> const
Operator::InitializeImplicationsMap() {
    return {{OperatorType::kEqual,
             {OperatorType::kEqual, OperatorType::kGreaterEqual, OperatorType::kLessEqual}},
            {OperatorType::kUnequal, {OperatorType::kUnequal}},
            {OperatorType::kGreater,
             {OperatorType::kGreater, OperatorType::kGreaterEqual, OperatorType::kUnequal}},
            {OperatorType::kLess,
             {OperatorType::kLess, OperatorType::kLessEqual, OperatorType::kUnequal}},
            {OperatorType::kGreaterEqual, {OperatorType::kGreaterEqual}},
            {OperatorType::kLessEqual, {OperatorType::kLessEqual}}};
}

std::unordered_map<OperatorType, std::vector<OperatorType>> const
Operator::InitializeTransitivesMap() {
    return {{OperatorType::kEqual, {OperatorType::kEqual}},
            {OperatorType::kUnequal, {OperatorType::kEqual}},
            {OperatorType::kGreater,
             {OperatorType::kGreater, OperatorType::kGreaterEqual, OperatorType::kEqual}},
            {OperatorType::kLess,
             {OperatorType::kLess, OperatorType::kLessEqual, OperatorType::kEqual}},
            {OperatorType::kGreaterEqual,
             {OperatorType::kGreater, OperatorType::kGreaterEqual, OperatorType::kEqual}},
            {OperatorType::kLessEqual,
             {OperatorType::kLess, OperatorType::kLessEqual, OperatorType::kEqual}}};
}

std::unordered_map<OperatorType, std::string> const Operator::InitializeShortStringMap() {
    return {{OperatorType::kEqual, "=="},        {OperatorType::kUnequal, "!="},
            {OperatorType::kGreater, ">"},       {OperatorType::kLess, "<"},
            {OperatorType::kGreaterEqual, ">="}, {OperatorType::kLessEqual, "<="}};
}

std::unordered_map<OperatorType, OperatorType> const Operator::kInverseMap =
        Operator::InitializeInverseMap();
std::unordered_map<OperatorType, OperatorType> const Operator::kSymmetricMap =
        Operator::InitializeSymmetricMap();
std::unordered_map<OperatorType, std::vector<OperatorType>> const Operator::kImplicationsMap =
        Operator::InitializeImplicationsMap();
std::unordered_map<OperatorType, std::vector<OperatorType>> const Operator::kTransitivesMap =
        Operator::InitializeTransitivesMap();
std::unordered_map<OperatorType, std::string> const Operator::kShortStringMap =
        Operator::InitializeShortStringMap();

size_t hash_value(model::Operator const& k) noexcept {
    return std::hash<model::Operator>()(k);
}

}  // namespace model
