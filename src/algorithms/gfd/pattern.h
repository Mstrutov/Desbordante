#pragma once
#include <Eigen/Dense>
#include <map>
#include <set>
#include <string>

class Pattern {
protected:
  std::set<int> indices;
  std::map<int, std::string> vertices;
  std::map<std::pair<int, int>, std::string> edges;

  Eigen::MatrixXi adjacencyMatrix;
  Eigen::MatrixXi degrees;
  int size;
  int vertexNum;

public:
  Pattern() = default;
  Pattern(const std::map<int, std::string> &,
          const std::map<std::pair<int, int>, std::string> &) noexcept(false);

  bool operator==(const Pattern &) const;
  bool operator!=(const Pattern &) const;

  std::set<int> getIndices() const { return this->indices; }
  std::map<int, std::string> getVertices() const { return this->vertices; }
  std::map<std::pair<int, int>, std::string> getEdges() const {
    return this->edges;
  }
  Eigen::MatrixXi getAdjacencyMatrix() const { return this->adjacencyMatrix; }
  Eigen::MatrixXi getVertexDegrees() const { return this->degrees; }
  int getSize() const { return this->size; }
  int getVertexNum() const { return this->vertexNum; }

  void print() const;
};
