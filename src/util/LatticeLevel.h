#pragma once

#include <map>
#include <vector>

#include "LatticeVertex.h"

class LatticeLevel{
private:
    unsigned int arity;
    std::map<boost::dynamic_bitset<>, std::unique_ptr<LatticeVertex>> vertices;
<<<<<<< HEAD

public:
    explicit LatticeLevel(unsigned int m_arity) : arity(m_arity) {}
    unsigned int getArity() const { return arity; }

=======
public:
    explicit LatticeLevel(unsigned int m_arity) : arity(m_arity) {}
    unsigned int getArity() const { return arity; }
>>>>>>> 81b41d2... TANE four-space indent fix

    std::map<boost::dynamic_bitset<>, std::unique_ptr<LatticeVertex>>& getVertices() { return vertices; }
    LatticeVertex const* getLatticeVertex(const boost::dynamic_bitset<>& columnIndices) const;
    void add(std::unique_ptr<LatticeVertex> vertex);

    //using vectors instead of lists because of .get()
    static void generateNextLevel(std::vector<std::unique_ptr<LatticeLevel>>& levels);
    static void clearLevelsBelow(std::vector<std::unique_ptr<LatticeLevel>>& levels, unsigned int arity);
};
