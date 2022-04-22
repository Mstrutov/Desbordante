#pragma once

#include <algorithm>
#include <string>

namespace util {

/* Class which represents string as vector of q-grams for calculating the cosine distance
 * between strings. Q-gram is a substring with length q. Q-gram vector is the vector of values,
 * which are the numbers of occurrences of each q-gram in the string. Cosine similarity is
 * calculated by taking the inner product of normalized q-gram vectors. Cosine distance is
 * calculated by subtracting cosine similarity from 1. For example: let q = 2, the string are:
 * "abc", "bcd". The q-gram vectors for the strings are: (1, 1, 0), (0, 1, 1). The first string
 * has 1 occurrence of "ab" and "bc" and 0 occurrences of "cd". Second string has 0 occurrences of
 * "ab" and 1 occurrence of "bc" and "cd". Cosine similarity between "abc" and "bcd" is equal to
 * (1*0 + 1*1 + 0*1) / (sqrt(1^2 + 1^2 + 0^2) * sqrt(1^2 + 1^2 + 0^2)) = 0.5.
 * Cosine distance between "abc" and "bcd" is equal to 1 - 0.5 = 0.5. */
class QGramVector {
private:
    double length_ = -1;
    std::unordered_map<std::string, unsigned> q_grams_;

    void CalculateLength();

public:
    explicit QGramVector(std::string_view const& string, unsigned q);

    double InnerProduct(QGramVector const& other) const;

    double GetLength() const {
        return length_;
    }

    double CosineSimilarity(QGramVector& other) const {
        return InnerProduct(other) / (GetLength() * other.GetLength());
    }

    double CosineDistance(QGramVector& other) const {
        return 1 - CosineSimilarity(other);
    }
};

}  // namespace util
