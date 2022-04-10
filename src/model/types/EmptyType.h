#pragma once

#include <cassert>

#include "Type.h"

namespace model {

class EmptyType : public virtual Type {
private:
    [[noreturn]] static inline void ThrowUnsupportedOperation() {
        throw std::logic_error("Meaningless operation");
    }

public:
    EmptyType() noexcept : Type(TypeId::kEmpty) {}

    [[nodiscard]]
    std::string ValueToString([[maybe_unused]] std::byte const* value) const override {
        return "";
    }

    CompareResult Compare([[maybe_unused]] std::byte const* l,
                          [[maybe_unused]] std::byte const* r) const override {
        ThrowUnsupportedOperation();
    }

    size_t Hash([[maybe_unused]] std::byte const* value) const override {
        ThrowUnsupportedOperation();
    }

    [[nodiscard]] size_t GetSize() const noexcept override {
        return 0;
    }

    void ValueFromStr([[maybe_unused]] std::byte* dest, std::string s) const override {
        if (!s.empty()) {
            throw std::invalid_argument("Cannot convert s to EmptyType value");
        }
    }

    [[nodiscard]] static std::byte* MakeValue() noexcept {
        return nullptr;
    }
};

}  // namespace model
