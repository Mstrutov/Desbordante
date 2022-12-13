#pragma once
#include <Eigen/Dense>
#include <map>
#include <set>
#include <string>

#include "pattern.h"

class Graph : public Pattern {
private:
  std::map<int, std::map<std::string, std::string>> attributes;

public:
  Graph() {};
  Graph(const std::map<int, std::string> &vertices_,
        const std::map<std::pair<int, int>, std::string> &edges_,
        const std::map<int, std::map<std::string, std::string>>
            &attributes_) noexcept(false)
      : Pattern(vertices_, edges_) {
    this->attributes = attributes_;
  }

  std::map<int, std::map<std::string, std::string>> getAttributes() const {
    return this->attributes;
  }

  void print() const;
};
