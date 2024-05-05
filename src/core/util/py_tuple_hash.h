#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>

namespace util {

class PyTupleHash {
    std::size_t res_ = 0x345678UL;
    std::size_t mult_ = 1000003UL;
    std::int64_t len_;

public:
    PyTupleHash(std::size_t len) noexcept : len_(len) {}

    [[nodiscard]] std::size_t GetResult() const noexcept {
        return res_;
    }

    void AppendHash(std::size_t hash) {
        --len_;
        res_ = (res_ ^ hash) * mult_;
        mult_ += 82520UL + len_ + len_;
    }

    template <typename T>
    void AddValue(T const& value) noexcept {
        AppendHash(std::hash<T>{}(value));
    }
};

}  // namespace util
