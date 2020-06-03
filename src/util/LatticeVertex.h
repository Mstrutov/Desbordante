// Strutovsky, 20.08.2019

#pragma once

#include <list>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "model/RelationalSchema.h"
#include "util/PositionListIndex.h"
//#include "../model/Column.h"
#include "../model/Vertical.h"


//enable_shared_from_this - if LV needs to create shared_ptr
class LatticeVertex{
private:
  Vertical vertical;
  std::shared_ptr<PositionListIndex> positionListIndex;
  dynamic_bitset<> rhsCandidates;
  bool isKeyCandidate = false;
  std::vector<std::shared_ptr<LatticeVertex>> parents;
  bool isInvalid = false;

  //list => vector due to usage of opertor[] in Tane
  //use pointer to vertical?
  // Using member initialization?
public:
  //TODO: no default initialization of PLI
  explicit LatticeVertex(Vertical&& _vertical) : vertical(_vertical), rhsCandidates(vertical.getSchema()->getNumColumns()) {}
  explicit LatticeVertex(Vertical& _vertical) : vertical(_vertical), rhsCandidates(vertical.getSchema()->getNumColumns()) {}

  std::vector<std::shared_ptr<LatticeVertex>>& getParents() { return parents; }
  //TODO: const - usually other Metanome classes use these outputs, so returning const isn't possible
  Vertical& getVertical() { return vertical; }
  boost::dynamic_bitset<>& getRhsCandidates() { return rhsCandidates; }

  void addRhsCandidates(std::vector<std::shared_ptr<Column>>&& candidates);

  //dynamic_bitset getBlockingPrefix();
  bool comesBeforeAndSharePrefixWith(LatticeVertex& that);
  bool getIsKeyCandidate() const { return isKeyCandidate; }
  void setKeyCandidate(bool m_isKeyCandidate) { isKeyCandidate = m_isKeyCandidate; }
  bool getIsInvalid() const { return isInvalid; }
  void setInvalid(bool m_isInvalid) { isInvalid = m_isInvalid; }

  //OK to store AND return ptr to PLI?
  std::shared_ptr<PositionListIndex> getPositionListIndex() { return positionListIndex; }
  void setPositionListIndex(shared_ptr<PositionListIndex> m_positionListIndex) { positionListIndex = m_positionListIndex; }

  //right analogy to compareTo?
  bool operator> (LatticeVertex& that);

  string toString();

  static bool comparator(const shared_ptr<LatticeVertex>& v1, const shared_ptr<LatticeVertex>& v2) { return *v2 > *v1; }
  friend std::ostream& operator<<(std::ostream& os, LatticeVertex& lv);
};
