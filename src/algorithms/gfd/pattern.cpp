#include <iostream>

#include "pattern.h"

Pattern::Pattern(
    const std::map<int, std::string> &vertices_,
    const std::map<std::pair<int, int>, std::string> &edges_) noexcept(false) {
  int max = 0;
  for (auto node : vertices_) {
    int index = node.first;
    this->indices.insert(index);
    max = max < index ? index : max;
  }
  if (max != this->indices.size() - 1) {
    throw std::out_of_range("index mismatch");
  }
  this->vertices = vertices_;
  this->edges = edges_;

  this->adjacencyMatrix =
      Eigen::MatrixXi::Zero(this->indices.size(), this->indices.size());
  for (auto edge_label : this->edges) {
    std::pair<int, int> edge = edge_label.first;
    this->adjacencyMatrix(edge.first, edge.second) = 1;
  }
  this->degrees = Eigen::MatrixXi::Zero(this->adjacencyMatrix.rows(), 2);
  for (int i = 0; i < this->adjacencyMatrix.rows(); ++i) {
    this->degrees(i, 0) = this->adjacencyMatrix.row(i).sum();
    this->degrees(i, 1) = this->adjacencyMatrix.col(i).sum();
  }
  this->vertexNum = this->vertices.size();
  this->size = vertexNum + this->edges.size();
}

bool Pattern::operator==(const Pattern &other) const {
  if (this->vertexNum != other.getVertexNum()) {
    return false;
  }
  if (this->size != other.getSize()) {
    return false;
  }
  Eigen::MatrixXi other_degrees = other.getVertexDegrees();
  Eigen::VectorXi check_out = Eigen::VectorXi::Zero(this->vertexNum);
  Eigen::VectorXi isomorphism = Eigen::VectorXi::Zero(this->vertexNum);
  for (int i = 0; i < this->vertexNum; ++i) {
    for (int j = 0; j < this->vertexNum; ++j) {
      if (this->degrees.row(i) == other_degrees.row(j) && !check_out(j)) {
        check_out(j) = 1;
        isomorphism(i) = j;
        break;
      }
    }
  }
  if (check_out.sum() != this->vertexNum) {
    return false;
  }
  Eigen::MatrixXi match =
      Eigen::MatrixXi::Zero(this->vertexNum, this->vertexNum);
  for (int i = 0; i < this->vertexNum; ++i) {
    match(i, isomorphism(i)) = 1;
  }
  Eigen::MatrixXi Q = this->adjacencyMatrix;
  Eigen::MatrixXi G = other.getAdjacencyMatrix();
  if ((match * (match * G).transpose()).transpose().cwiseProduct(Q) != Q) {
    return false;
  }
  Eigen::VectorXi mask = Eigen::VectorXi::Zero(G.rows());
  for (int i = 0; i < G.rows(); ++i) {
    mask(i) = i;
  }
  for (const auto &edge_label : this->edges) {
    std::pair<int, int> edge = edge_label.first;
    Eigen::VectorXi index1 = Eigen::VectorXi::Zero(Q.rows());
    index1(edge.first) = 1;
    int i = mask.cwiseProduct((index1.transpose() * match).transpose()).sum();
    Eigen::VectorXi index2 = Eigen::VectorXi::Zero(Q.rows());
    index2(edge.second) = 1;
    int j = mask.cwiseProduct((index2.transpose() * match).transpose()).sum();
    if (edge_label.second != other.getEdges().at(std::pair<int, int>(i, j))) {
      return false;
    }
  }
  return true;
}

bool Pattern::operator!=(const Pattern &other) const {
  return !(*this == other);
}

void Pattern::print() const {
  std::cout << "Adjacency matrix:" << std::endl
            << this->adjacencyMatrix << std::endl;
  std::cout << "Vertex labels:" << std::endl;
  for (const auto &node : this->vertices) {
    std::cout << node.first << ":" << node.second << "; ";
  }
  std::cout << std::endl << "Edge labels:" << std::endl;
  for (const auto &edge_label : this->edges) {
    auto edge = edge_label.first;
    std::cout << "(" << edge.first << "," << edge.second
              << "):" << edge_label.second << "; ";
  }
}
