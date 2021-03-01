#include <iostream>

#pragma once

// ~ Java.Random
class CustomRandom {
private:
    long long seed = 12345;
    long long addend = 0xBLL;
    long long multiplier = 0x5DEECE66DLL;
    long long mask = (1LL << 48) - 1;

    const double DOUBLE_UNIT = 1.0 / (1LL << 53);

    long long abs(long long number) {
        if (number < 0) {
            return -number;
        } else {
            return number;
        }
    }
public:
    explicit CustomRandom(long long seed = 12345) : seed(initialScramble(seed)) {
        //std::cout << "Initializing with seed = " << seed << '\n';
    }

    static long long initialScramble(long long seed) { return (seed ^ 25214903917LL) & 281474976710655LL; }

    long long next(int bits) {
        seed = (seed * multiplier + addend) & mask;
        return (int)(static_cast<unsigned long long>(seed) >> (48 - bits));
    }

    long long nextLL() {
        long long res = ((long long)(next(32)) << 32) + next(32);
        return res; 
    }

    long long nextULL() {
        return abs(nextLL());
    }

    int nextInt() {
        int res = next(32);
        // std::cout << "nextInt generated " <<  res << '\n';
        return res;
    }

    int nextInt(int upperBound) {
        if (upperBound <= 0) {
            throw std::invalid_argument("UpperBound must be >= 0");
        }

        int res = next(31);
        int m = upperBound - 1;
        if ((upperBound & m) == 0)
            res = (int) ((upperBound * (long long)res) >> 31);
        else {
            for (int u = res;
                 u - (res = u % upperBound) + m < 0;
                 u = next(31))
                ;
        }
        // std::cout << "bounded nextInt generated " <<  res << '\n';
        return res;
    }

    double nextDouble() { 
        double res = (((long long)(next(26)) << 27) + next(27)) * DOUBLE_UNIT;
        // std::cout << "nextDouble generated " << res << '\n';
        return res; 
    }
};
