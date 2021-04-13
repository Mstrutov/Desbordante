//
// Created by alexandrsmirn
//

#pragma once

#include "Vertical.h"
#include "PositionListIndex.h"

enum class NodeCategory {
    unvisited,
    dependency,
    minimalDependency,
    candidateMinimalDependency,
    nonDependency,
    maximalNonDependency,
    candidateMaximalNonDependency
};

class LatticeNode : Vertical{
private:
    //Vertical vertical; -- наследование или агрегирование?
    NodeCategory category;
    shared_ptr<PositionListIndex> partition;

public:
    LatticeNode(Column const& column);

    LatticeNode() = delete;
    LatticeNode(LatticeNode const& other) = delete;
    LatticeNode& operator=(LatticeNode const& other) = delete;
    LatticeNode(LatticeNode&& other) = delete;
    LatticeNode& operator=(LatticeNode&& other) = delete;
    virtual ~LatticeNode() = default; //TODO: подумать про методы и особенно про деструктор


    bool isVisited() { return category != NodeCategory::unvisited; }
    bool isCandidate();
    bool isMinimalDependency();
    bool isMaximalNonDependency();
    //TODO: еще какие-нибудь методы
};
