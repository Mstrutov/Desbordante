#pragma once

#include "Type.h"

namespace model {

class StringType : public Type {
public:
    /* type_id parameter is temporary for BigIntType */
    explicit StringType(TypeId type_id = TypeId::kString) noexcept : Type(type_id) {}

    void Free(std::byte const* value) const noexcept override {
        Destruct(value);
        Type::Free(value);
    }

    [[nodiscard]] std::string ValueToString(std::byte const* value) const override {
        return GetValue<String>(value);
    }

    [[nodiscard]] CompareResult Compare(std::byte const* l, std::byte const* r) const override {
        auto const& l_val = GetValue<String>(l);
        auto const& r_val = GetValue<String>(r);

        int const res = l_val.compare(r_val);
        if (res == 0) {
            return CompareResult::kEqual;
        }
        if (res < 0) {
            return CompareResult::kLess;
        }
        return CompareResult::kGreater;
    }

    [[nodiscard]] size_t Hash(std::byte const* value) const override {
        return std::hash<String>()(GetValue<String>(value));
    }

    [[nodiscard]] size_t GetSize() const override {
        return sizeof(String);
    }

    void ValueFromStr(std::byte* dest, std::string s) const override {
        new (dest) String(std::move(s));
    }

    [[nodiscard]] std::byte* MakeValue(String v) const {
        std::byte* buf = Allocate();
        new (buf) String(std::move(v));
        return buf;
    }

    static void Destruct(std::byte const* v) {
        reinterpret_cast<String const*>(v)->~String();
    }
};

}  // namespace model
