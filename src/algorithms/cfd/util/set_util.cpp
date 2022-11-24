#include <algorithm>
#include <numeric>

// see ../algorithms/cfd/LICENSE

#include "algorithms/cfd/util/set_util.h"

[[maybe_unused]] SubsetIterator::SubsetIterator(int size, int subSize)
: seed((1 << subSize) - 1), max_seed(1 << size), subs_number(BinomialCoeff(size, subSize)) {

}

[[maybe_unused]] std::bitset<64> SubsetIterator::Next() {
    int x = seed;
    int u = x & (-x);
    int v = x + u;
    seed = v + (((v ^ x) / u) >> 2);
    return std::bitset<64>(x);
}

[[maybe_unused]] long int SubsetIterator::GetSubsNumber() {
    return subs_number;
}

[[maybe_unused]] std::vector<int> Range(int min, int max, int step) {
    std::vector<int> res;
    res.reserve((max - min) / step);
    for (int i = min; i != max; i += step) {
        res.push_back(i);
    }
    return res;
}

// Возвращает список вида 0 1 2 ... размера max
[[maybe_unused]] std::vector<int> Iota(int max) {
    std::vector<int> iotas(max);
    std::iota(iotas.begin(), iotas.end(), 0);
    return iotas;
}

[[maybe_unused]] int BinomialCoeff(int n, int k)
{
    int res = 1;
 
    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;
 
    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }
 
    return res;
}

[[maybe_unused]] std::vector<std::vector<int>> CartesianProduct(const std::vector<int>& sizes) {
    std::vector<int> counters(sizes.size());
    std::vector<std::vector<int>> cart_prod;
    cart_prod.reserve(product(sizes));

    while (counters[0] < sizes[0]) {
        std::vector<int> row(sizes.size());
        for (unsigned i = 0; i < sizes.size(); i++) {
            row[i] = counters[i];
        }
        int inc_ix = counters.size() - 1;
        counters[inc_ix]++;
        while (inc_ix > 0 && counters[inc_ix] == sizes[inc_ix]) {
            counters[inc_ix] = 0;
            counters[--inc_ix]++;
        }
        cart_prod.push_back(row);
    }
    return cart_prod;
}

[[maybe_unused]] void SubsetsLengthK(const int size, const int k, std::vector<std::bitset<32>>& subs) {
    const int max_b = pow(2, size);
    int x = pow(2, k) - 1;
    while (x < max_b) {
        subs.push_back(std::bitset<32>(x));
        int u = x & (-x);
        int v = x + u;
        x = v + (((v ^ x) / u) >> 2);
    }
}

[[maybe_unused]] std::vector<std::bitset<32>> AllSubsets(const int size) {
    std::vector<std::bitset<32>> subs;
    subs.reserve(pow(2, size)-2);
    for (int i = 1; i < size; i++) {
        SubsetsLengthK(size, i, subs);
    }
    return subs;
}

[[maybe_unused]] std::vector<std::bitset<32>> AllSubsetsIncl(const int size) {
    std::vector<std::bitset<32>> subs;
    subs.reserve(pow(2, size)-2);
    for (int i = 1; i <= size; i++) {
        SubsetsLengthK(size, i, subs);
    }
    return subs;
}
