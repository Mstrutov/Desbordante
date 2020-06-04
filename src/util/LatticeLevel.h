// Strutovsky, 21.08.2019

#pragma once

#include <map>
#include <vector>

#include "LatticeVertex.h"

class LatticeLevel{
private:
  int arity;
  //TODO: store by value only if LatticeVertex exists only as a member of Level
  //vertices are used once to get some data from LatticeLevel
  //but in main algorithm vertices are firstly created, then added to a level =>
  //perhaps, use unique + std::move
  std::map<boost::dynamic_bitset<>, std::shared_ptr<LatticeVertex>> vertices;
public:
  explicit LatticeLevel(int m_arity) { arity = m_arity; }
  int getArity() const { return arity; }

  std::map<boost::dynamic_bitset<>, std::shared_ptr<LatticeVertex>>& getVertices() { return vertices; }
  std::shared_ptr<LatticeVertex> getLatticeVertex(const boost::dynamic_bitset<>& columnIndices);
  void add(std::shared_ptr<LatticeVertex> vertex);

  //using vectors instead of lists because of .get()
  static void generateNextLevel(std::vector<std::shared_ptr<LatticeLevel>>& levels);
  static void clearLevelsBelow(std::vector<std::shared_ptr<LatticeLevel>>& levels, int arity);
};
