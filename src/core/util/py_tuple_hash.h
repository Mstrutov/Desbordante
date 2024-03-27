#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>

namespace util {

template <typename T>
class PyTupleHash {
    std::size_t res_ = 0x345678UL;
    std::size_t mult_ = 1000003UL;
    std::hash<T> const hasher_{};
    std::int64_t len_;

public:
    PyTupleHash(std::size_t len) noexcept : len_(len) {}

    [[nodiscard]] std::size_t GetResult() const noexcept {
        return res_;
    }

    void AddValue(T const& value) noexcept {
        --len_;
        std::size_t hash = hasher_(value);
        res_ = (res_ ^ hash) * mult_;
        mult_ += 82520UL + len_ + len_;
    }
};

}  // namespace util
