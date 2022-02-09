#pragma once

#include <map>
#include <vector>

#include "LatticeVertex.h"

namespace util {

class LatticeLevel {
private:
    unsigned int arity_;
    std::map<boost::dynamic_bitset<>, std::unique_ptr<LatticeVertex>> vertices_;

public:
    explicit LatticeLevel(unsigned int m_arity) : arity_(m_arity) {}
    unsigned int GetArity() const { return arity_; }

    std::map<boost::dynamic_bitset<>, std::unique_ptr<LatticeVertex>>& GetVertices() {
        return vertices_;
    }
    LatticeVertex const* GetLatticeVertex(const boost::dynamic_bitset<>& column_indices) const;
    void Add(std::unique_ptr<LatticeVertex> vertex);

    //using vectors instead of lists because of .get()
    static void GenerateNextLevel(std::vector<std::unique_ptr<LatticeLevel>>& levels);
    static void ClearLevelsBelow(std::vector<std::unique_ptr<LatticeLevel>>& levels,
                                 unsigned int arity);
};

} // namespace util

