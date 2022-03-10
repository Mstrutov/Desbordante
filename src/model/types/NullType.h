#pragma once

#include <atomic>

#include "Type.h"

namespace model {

class NullType : public virtual Type {
private:
    bool is_null_eq_null_;

public:
    explicit NullType(bool is_null_eq_null) noexcept
        : Type(TypeId::kNull), is_null_eq_null_(is_null_eq_null) {}

    [[nodiscard]] bool IsNullEqNull() const noexcept {
        return is_null_eq_null_;
    }

    [[nodiscard]]
    std::string ValueToString([[maybe_unused]] std::byte const* value) const override {
        return Null::kValue.data();
    }

    [[nodiscard]] CompareResult Compare([[maybe_unused]] std::byte const* l,
                                        [[maybe_unused]] std::byte const* r) const override {
        if (is_null_eq_null_) {
            return CompareResult::kEqual;
        }
        return CompareResult::kNotEqual;
    }

    [[nodiscard]] size_t Hash([[maybe_unused]] std::byte const* value) const override {
        static std::atomic_uint hash = 0;
        if (is_null_eq_null_) {
            return 1;
        }
        return hash.fetch_add(1, std::memory_order_relaxed);
    }

    [[nodiscard]] size_t GetSize() const override {
        return 0;
    }

    void ValueFromStr([[maybe_unused]] std::byte* dest, std::string s) const override {
        if (s != Null::kValue) {
            throw std::invalid_argument("Cannot convert s to NullType value");
        }
    }

    [[nodiscard]] static std::byte* MakeValue() noexcept {
        return nullptr;
    }
};

}  // namespace model
