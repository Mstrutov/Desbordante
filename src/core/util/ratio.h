#pragma once

namespace util {
template <typename Type>
struct Ratio {
    Type numerator;
    Type denominator;

    constexpr Ratio(Type numerator, Type denominator)
        : numerator(numerator), denominator(denominator) {}

    Ratio& operator *=(Ratio other) {
        numerator *= other.numerator;
        denominator *= other.denominator;
        return *this;
    }
};
}  // namespace util
