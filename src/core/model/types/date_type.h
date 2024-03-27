#pragma once
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>

#include <boost/date_time/gregorian/gregorian.hpp>

#include "imetrizable_type.h"
#include "type.h"

namespace model {

class DateType : public IMetrizableType {
public:
    using Delta = boost::gregorian::date_duration;

    DateType() noexcept : IMetrizableType(TypeId::kDate) {}

    void Free(std::byte const* value) const noexcept override {
        Destruct(value);
        Type::Free(value);
    }

    [[nodiscard]] std::string ValueToString(std::byte const* value) const override {
        return boost::gregorian::to_iso_extended_string(GetValue<Date>(value));
    }

    [[nodiscard]] size_t Hash(std::byte const* value) const override {
        return GetValue<Date>(value).julian_day();
    }

    void ValueFromStr(std::byte* dest, std::string s) const override {
        new (dest) Date(boost::gregorian::from_simple_string(s));
    }

    [[nodiscard]] std::unique_ptr<Type> CloneType() const override {
        return std::make_unique<DateType>();
    }

    [[nodiscard]] CompareResult Compare(std::byte const* l, std::byte const* r) const override {
        auto const& l_val = GetValue<Date>(l);
        auto const& r_val = GetValue<Date>(r);
        return Compare(l_val, r_val);
    }

    static CompareResult Compare(Date const& l_val, Date const& r_val) {
        if (l_val == r_val) {
            return CompareResult::kEqual;
        }
        if (l_val < r_val) {
            return CompareResult::kLess;
        }
        return CompareResult::kGreater;
    }

    [[nodiscard]] size_t GetSize() const override {
        return sizeof(Date);
    }

    [[nodiscard]] std::byte* MakeValue(Date const& date = Date()) const {
        auto* buf = Allocate();
        new (buf) Date(date);
        return buf;
    }

    std::byte* AddDelta(std::byte const* date, Delta const& delta, std::byte* res) const {
        GetValue<Date>(res) = GetValue<Date>(date) + delta;
        return res;
    }

    std::byte* AddDeltaFromInt(std::byte const* date, std::byte const* delta,
                               std::byte* res) const {
        return AddDelta(date, Delta(GetValue<Int>(delta)), res);
    }

    std::byte* AddDelta(std::byte const* date, Delta const& delta) const {
        std::byte* res = MakeValue();
        return AddDelta(date, delta, res);
    }

    std::byte* SubDelta(std::byte const* date, Delta const& delta, std::byte* res) const {
        GetValue<Date>(res) = GetValue<Date>(date) - delta;
        return res;
    }

    std::byte* SubDeltaFromInt(std::byte const* date, std::byte const* delta,
                               std::byte* res) const {
        return SubDelta(date, Delta(GetValue<Int>(delta)), res);
    }

    std::byte* SubDelta(std::byte const* date, Delta const& delta) const {
        std::byte* res = MakeValue();
        return SubDelta(date, delta, res);
    }

    Delta SubDate(std::byte const* left_date, std::byte const* right_date) const {
        Delta res = GetValue<Date>(left_date) - GetValue<Date>(right_date);
        return res;
    }

    std::byte* SubDate(std::byte const* left_date, std::byte const* right_date,
                       std::byte* res) const {
        GetValue<Int>(res) = SubDate(left_date, right_date).days();
        return res;
    }

    double Dist(std::byte const* l, std::byte const* r) const override {
        Delta d = SubDate(l, r);
        return std::abs(d.days());
    }

    static void Destruct(std::byte const* value) {
        GetValue<Date>(value).~Date();
    }

    [[nodiscard]] auto GetDeleater() const {
        return [this](std::byte const* v) { Free(v); };
    }

    [[nodiscard]] Destructor GetDestructor() const override {
        return Destruct;
    }
};

using DateTypeDeleter = decltype(std::declval<DateType>().GetDeleater());

}  // namespace model
