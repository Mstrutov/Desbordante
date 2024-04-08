#pragma once

#include "imetrizable_type.h"
#include "type.h"
#include "util/levenshtein_distance.h"

namespace model {

class StringType : public IMetrizableType {
public:
    /* type_id parameter is temporary for BigIntType */
    explicit StringType(TypeId type_id = TypeId::kString) noexcept : IMetrizableType(type_id) {}

    void Free(std::byte const* value) const noexcept override {
        Destruct(value);
        Type::Free(value);
    }

    [[nodiscard]] std::string ValueToString(std::byte const* value) const override {
        return GetValue<String>(value);
    }

    [[nodiscard]] std::byte* Clone(std::byte const* value) const override {
        return MakeValue(ValueToString(value));
    }

    [[nodiscard]] std::unique_ptr<Type> CloneType() const override {
        return std::make_unique<StringType>();
    }

    [[nodiscard]] CompareResult Compare(std::byte const* l, std::byte const* r) const override {
        auto const& l_val = GetValue<String>(l);
        auto const& r_val = GetValue<String>(r);

        return Compare(l_val, r_val);
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

    [[nodiscard]] std::byte* MakeValue(String v = "") const {
        std::byte* buf = Allocate();
        new (buf) String(std::move(v));
        return buf;
    }

    static CompareResult Compare(String const& l_val, String const& r_val) {
        int const res = l_val.compare(r_val);
        if (res == 0) {
            return CompareResult::kEqual;
        }
        if (res < 0) {
            return CompareResult::kLess;
        }
        return CompareResult::kGreater;
    }

    std::byte* Concat(std::byte* l, std::byte* r, std::byte* result) const {
        String const& l_val = GetValue<String>(l);
        String const& r_val = GetValue<String>(r);

        GetValue<String>(result) = l_val + r_val;
        return result;
    }

    std::byte* Concat(std::byte* l, std::byte* r) const {
        std::byte* result = MakeValue();
        return Concat(l, r, result);
    }

    auto GetDeleter() const {
        return [this](std::byte const* v) { Free(v); };
    }

    /* Returns Levenshtein distance between l and r */
    double Dist(std::byte const* l, std::byte const* r) const override {
        return util::LevenshteinDistance(GetValue<String>(l), GetValue<String>(r));
    }

    static void Destruct(std::byte const* v) {
        reinterpret_cast<String const*>(v)->~String();
    }

    Destructor GetDestructor() const override {
        return Destruct;
    }
};

using StringTypeDeleter = decltype(std::declval<StringType>().GetDeleter());

}  // namespace model
