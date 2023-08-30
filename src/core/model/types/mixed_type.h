#pragma once

#include <cassert>
#include <memory>
#include <unordered_map>

#include "big_int_type.h"
#include "create_type.h"
#include "date_type.h"
#include "double_type.h"
#include "empty_type.h"
#include "int_type.h"
#include "null_type.h"
#include "numeric_type.h"
#include "string_type.h"
#include "type.h"
#include "undefined_type.h"

namespace model {

/* Holds value of any type of TypeId.
 * First kTypeIdSize bytes of std::byte* value representation should be a TypeId of value,
 * all other bytes represent the value itself
 */
class MixedType final : public Type {
private:
    bool is_null_eq_null_;

    [[nodiscard]] static size_t GetTypeIdSizeWithPadding(TypeId type_id) noexcept {
        size_t alignment = GetAlignment(type_id);
        assert(kTypeIdSize <= alignment);
        return alignment;
    }

    static constexpr size_t kTypeIdSize = sizeof(TypeId::_integral);

public:
    explicit MixedType(bool is_null_eq_null) noexcept
        : Type(TypeId::kMixed), is_null_eq_null_(is_null_eq_null) {}

    [[nodiscard]] std::string ValueToString(std::byte const* value) const final {
        std::unique_ptr<Type> type = RetrieveType(value);
        return type->ValueToString(RetrieveValue(value));
    }

    [[nodiscard]] std::unique_ptr<Type> CloneType() const override {
        return std::make_unique<MixedType>(is_null_eq_null_);
    }

    [[nodiscard]] CompareResult Compare(std::byte const* l, std::byte const* r) const final {
        TypeId l_type_id = RetrieveTypeId(l);
        TypeId r_type_id = RetrieveTypeId(r);

        if (l_type_id != r_type_id) {
            throw std::invalid_argument("Cannot compare values of different types");
        }

        std::unique_ptr<Type> type = RetrieveType(l);
        return type->Compare(RetrieveValue(l), RetrieveValue(r));
    }

    [[nodiscard]] CompareResult CompareAsStrings(std::byte const* l, std::byte const* r) const {
        String l_str = this->ValueToString(l);
        String r_str = this->ValueToString(r);

        return StringType::Compare(l_str, r_str);
    }

    void Free(std::byte const* value) const noexcept override {
        TypeId const type_id = RetrieveTypeId(value);
        if (type_id == +TypeId::kString || type_id == +TypeId::kBigInt) {
            StringType::Destruct(RetrieveValue(value));
        }
        if (type_id == +TypeId::kDate) {
            DateType::Destruct(RetrieveValue(value));
        }
        Type::Free(value);
    }

    [[nodiscard]] size_t Hash(std::byte const* value) const final {
        std::unique_ptr<Type> type = RetrieveType(value);
        return type->Hash(RetrieveValue(value));
    }

    [[nodiscard]] size_t GetSize() const final {
        throw std::logic_error("Mixed type does not have a fixed size");
    }

    // It's correct, but not optimal, need to be rewrited later with other virtual
    // Clone(std::byte const* value, std::byte const* new_value)
    [[nodiscard]] std::byte* Clone(std::byte const* value) const override {
        std::unique_ptr<Type> type = RetrieveType(value);
        size_t size = GetMixedValueSize(type.get());
        auto* new_value = new std::byte[size];
        ValueFromStr(new_value, ValueToString(value), type.get());
        return new_value;
    }

    /* Note: first kTypeIdSize bytes of dest should contain type_id to cast to */
    void ValueFromStr(std::byte* dest, std::string s) const final {
        std::unique_ptr<Type> type = RetrieveType(dest);
        type->ValueFromStr(RetrieveValue(dest), s);
    }

    [[nodiscard]] bool IsNullEqNull() const noexcept {
        return is_null_eq_null_;
    }

    template <typename T>
    [[nodiscard]] std::byte* MakeValue(T literal, Type const* type) const {
        std::byte* buf = AllocateMixed(type);
        new (RetrieveValue(buf)) T(std::move(literal));
        return buf;
    }

    [[nodiscard]] std::unique_ptr<Type> RetrieveType(std::byte const* value) const {
        TypeId type_id = RetrieveTypeId(value);
        return CreateType(type_id, is_null_eq_null_);
    }

    [[nodiscard]] static size_t GetAlignment(TypeId type_id) noexcept {
        // Could possibly be implemented as a virtual method, however, it's much more incovinient
        // since virtual method requires an object of type `Type` to be called which is not always
        // available. Moreover, alignment of a type is used only when implementing MixedType
        // or working with it, so GetAlignment as a general method for any type seems quite useless.
        static std::unordered_map<TypeId, size_t> const type_to_alignment = {
                {TypeId::kNull, alignof(Null)},     {TypeId::kEmpty, alignof(Empty)},
                {TypeId::kInt, alignof(Int)},       {TypeId::kDouble, alignof(Double)},
                {TypeId::kString, alignof(String)}, {TypeId::kBigInt, alignof(BigInt)}};
        assert(type_to_alignment.find(type_id) != type_to_alignment.end());
        return type_to_alignment.at(type_id);
    }

    static void ValueFromStr(std::byte* dest, std::string s, Type const* type) {
        type->ValueFromStr(SetTypeId(dest, type->GetTypeId()), std::move(s));
    }

    [[nodiscard]] static std::byte* AllocateMixed(Type const* type) {
        auto* buf = new std::byte[GetMixedValueSize(type)];
        RetrieveTypeId(buf) = type->GetTypeId();
        return buf;
    }

    [[nodiscard]] static TypeId const& RetrieveTypeId(std::byte const* value) noexcept {
        return GetValue<TypeId>(value);
    }

    [[nodiscard]] static TypeId& RetrieveTypeId(std::byte* value) noexcept {
        return GetValue<TypeId>(value);
    }

    [[nodiscard]] static std::byte* RetrieveValue(std::byte* value_with_type) noexcept {
        return value_with_type + GetTypeIdSizeWithPadding(RetrieveTypeId(value_with_type));
    }

    [[nodiscard]] static std::byte const* RetrieveValue(std::byte const* value_with_type) noexcept {
        return value_with_type + GetTypeIdSizeWithPadding(RetrieveTypeId(value_with_type));
    }

    static std::byte* SetTypeId(std::byte* dest, TypeId const type_id) noexcept {
        RetrieveTypeId(dest) = type_id;
        return RetrieveValue(dest);
    }

    [[nodiscard]] static size_t GetMixedValueSize(Type const* type) noexcept {
        return GetTypeIdSizeWithPadding(type->GetTypeId()) + type->GetSize();
    }
};

}  // namespace model
