#pragma once

#include "EmptyType.h"
#include "NullType.h"

namespace model {

class UndefinedType final : public EmptyType, public NullType {
private:
    [[noreturn]] static inline void ThrowUnsupportedOperation() {
        throw std::logic_error("Meaningless operation. Use EmptyType or NullType methods directly");
    }

public:
    explicit UndefinedType(bool is_null_equal_null) noexcept
        : Type(TypeId::kUndefined), NullType(is_null_equal_null) {}

    std::string ValueToString([[maybe_unused]] std::byte const* value) const final {
        ThrowUnsupportedOperation();
    }

    CompareResult Compare([[maybe_unused]] std::byte const* l,
                          [[maybe_unused]] std::byte const* r) const final {
        ThrowUnsupportedOperation();
    }

    size_t Hash([[maybe_unused]] std::byte const* value) const final {
        ThrowUnsupportedOperation();
    }

    [[nodiscard]] size_t GetSize() const noexcept final {
        return 0;
    }

    void ValueFromStr([[maybe_unused]] std::byte* dest,
                      [[maybe_unused]] std::string s) const final {
        ThrowUnsupportedOperation();
    }
};

}  // namespace model
