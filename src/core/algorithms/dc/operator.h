#pragma once

#include <boost/functional/hash.hpp>

#include "type.h"

namespace model {

enum class OperatorType { kEqual, kUnequal, kGreater, kLess, kGreaterEqual, kLessEqual };

constexpr std::array<OperatorType, 6> kAllOperatorTypes = {
        OperatorType::kEqual, OperatorType::kUnequal,      OperatorType::kGreater,
        OperatorType::kLess,  OperatorType::kGreaterEqual, OperatorType::kLessEqual};

class Operator {
private:
    OperatorType op_;
    static std::unordered_map<OperatorType, OperatorType> const kInverseMap;
    static std::unordered_map<OperatorType, OperatorType> const kSymmetricMap;
    static std::unordered_map<OperatorType, std::vector<OperatorType>> const kImplicationsMap;
    static std::unordered_map<OperatorType, std::vector<OperatorType>> const kTransitivesMap;
    static std::unordered_map<OperatorType, std::string> const kShortStringMap;

public:
    Operator(OperatorType type) : op_(type) {}

    bool operator==(Operator const& rhs) const {
        return op_ == rhs.op_;
    }

    bool Eval(std::byte const* v1, std::byte const* v2, Type const& type) const;

    // 'a op b' <=> !'a op.inverse b'
    Operator GetInverse() const {
        return Operator(kInverseMap.at(op_));
    }

    // 'a op b' <=> 'b op.symmetric a'
    Operator GetSymmetric() const {
        return Operator(kSymmetricMap.at(op_));
    }

    // If 'a op b', then 'a op.implications[i] b'
    std::vector<Operator> GetImplications() const {
        std::vector<Operator> implications;
        for (auto type : kImplicationsMap.at(op_)) {
            implications.push_back(Operator(type));
        }
        return implications;
    }

    // If 'a op b' and 'b op.transitives[i] c', then 'a op c'
    std::vector<Operator> GetTransitives() const {
        std::vector<Operator> transitives;
        for (auto type : kTransitivesMap.at(op_)) {
            transitives.push_back(Operator(type));
        }
        return transitives;
    }

    std::string ToString() const {
        return kShortStringMap.at(op_);
    }

    OperatorType GetType() const {
        return op_;
    }

    static std::unordered_map<OperatorType, OperatorType> const InitializeInverseMap();

    static std::unordered_map<OperatorType, OperatorType> const InitializeSymmetricMap();

    static std::unordered_map<OperatorType, std::vector<OperatorType>> const
    InitializeImplicationsMap();

    static std::unordered_map<OperatorType, std::vector<OperatorType>> const
    InitializeTransitivesMap();

    static std::unordered_map<OperatorType, std::string> const InitializeShortStringMap();
};

// NOLINTBEGIN(readability-identifier-naming)
size_t hash_value(model::Operator const& k) noexcept;
// NOLINTEND(readability-identifier-naming)

}  // namespace model

namespace std {
template <>
struct hash<model::Operator> {
    size_t operator()(model::Operator const& k) const noexcept {
        return static_cast<size_t>(k.GetType());
    }
};
};  // namespace std
