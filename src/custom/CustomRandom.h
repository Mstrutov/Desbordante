#include <iostream>

#pragma once

// ~ Java.Random
class CustomRandom {
private:
    long long seed_ = 12345;
    long long addend_ = 0xBLL;
    long long multiplier_ = 0x5DEECE66DLL;
    long long mask_ = (1LL << 48) - 1;

    const double kDoubleUnit = 1.0 / (1LL << 53);

    long long Abs(long long number) {
        if (number < 0) {
            return -number;
        } else {
            return number;
        }
    }
public:
    explicit CustomRandom(long long seed = 12345) : seed_(InitialScramble(seed)) {
        //std::cout << "Initializing with seed = " << seed << '\n';
    }

    static long long InitialScramble(long long seed) {
        return (seed ^ 25214903917LL) & 281474976710655LL;
    }

    long long Next(int bits) {
        seed_ = (seed_ * multiplier_ + addend_) & mask_;
        return (int)(static_cast<unsigned long long>(seed_) >> (48 - bits));
    }

    long long NextLL() {
        long long res = ((long long)(Next(32)) << 32) + Next(32);
        return res;
    }

    long long NextULL() {
        return Abs(NextLL());
    }

    int NextInt() {
        int res = Next(32);
        // std::cout << "NextInt generated " <<  res << '\n';
        return res;
    }

    int NextInt(int upper_bound) {
        if (upper_bound <= 0) {
            throw std::invalid_argument("UpperBound must be >= 0");
        }

        int res = Next(31);
        int m = upper_bound - 1;
        if ((upper_bound & m) == 0)
            res = (int)((upper_bound * (long long)res) >> 31);
        else {
            for (int u = res; u - (res = u % upper_bound) + m < 0; u = Next(31))
                ;
        }
        // std::cout << "bounded NextInt generated " <<  res << '\n';
        return res;
    }

    double NextDouble() {
        double res = (((long long)(Next(26)) << 27) + Next(27)) * kDoubleUnit;
        // std::cout << "NextDouble generated " << res << '\n';
        return res;
    }
};
