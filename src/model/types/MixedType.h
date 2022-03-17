#pragma once

#include <cassert>
#include <memory>

#include "BigIntType.h"
#include "CreateType.h"
#include "DoubleType.h"
#include "EmptyType.h"
#include "IntType.h"
#include "NullType.h"
#include "NumericType.h"
#include "StringType.h"
#include "Type.h"
#include "UndefinedType.h"

namespace model {

/* Holds value of any type of TypeId.
 * First kTypeIdSize bytes of std::byte* value representation should be a TypeId of value,
 * all other bytes represent the value itself
 */
class MixedType final : public Type {
private:
    bool is_null_eq_null_;

public:
    static constexpr size_t kTypeIdSize = sizeof(TypeId::_integral);

    explicit MixedType(bool is_null_eq_null) noexcept
        : Type(TypeId::kMixed), is_null_eq_null_(is_null_eq_null) {}

    [[nodiscard]] std::string ValueToString(std::byte const* value) const final {
        std::unique_ptr<Type> type = RetrieveType(value);
        return type->ValueToString(RetrieveValue(value));
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

    void Free(std::byte const* value) const noexcept override {
        TypeId const type_id = RetrieveTypeId(value);
        if (type_id == +TypeId::kString || type_id == +TypeId::kBigInt) {
            StringType::Destruct(RetrieveValue(value));
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

    [[nodiscard]] std::byte* Clone(std::byte const* value) const override {
        std::unique_ptr<Type> type = RetrieveType(value);
        size_t size = GetMixedValueSize(type.get());
        auto* new_value = new std::byte[size];
        std::memcpy(new_value, value, size);
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

    static void ValueFromStr(std::byte* dest, std::string s, Type const* type) {
        type->ValueFromStr(SetTypeId(dest, type->GetTypeId()), std::move(s));
    }

    [[nodiscard]] static std::byte* AllocateMixed(Type const* type) {
        auto* buf = new std::byte[GetMixedValueSize(type)];
        *reinterpret_cast<TypeId*>(buf) = type->GetTypeId();
        return buf;
    }

    [[nodiscard]] static TypeId RetrieveTypeId(std::byte const* value) {
        TypeId type_id = *reinterpret_cast<TypeId const*>(value);
        return type_id;
    }

    [[nodiscard]] static std::byte* RetrieveValue(std::byte* value_with_type) noexcept {
        return value_with_type + kTypeIdSize;
    }

    [[nodiscard]] static std::byte const* RetrieveValue(std::byte const* value_with_type) noexcept {
        return value_with_type + kTypeIdSize;
    }

    static std::byte* SetTypeId(std::byte* dest, TypeId const type_id) {
        *reinterpret_cast<TypeId*>(dest) = type_id;
        return dest + kTypeIdSize;
    }

    [[nodiscard]] static size_t GetMixedValueSize(Type const* type) noexcept {
        return kTypeIdSize + type->GetSize();
    }
};

}  // namespace model
