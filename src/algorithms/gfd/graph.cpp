#include <iostream>

#include "graph.h"

void Graph::print() const {
  Pattern::print();
  std::cout << std::endl << "Attributes:";
  for (const auto &attr : this->attributes) {
    std::cout << std::endl << attr.first << " -> ";
    for (const auto &value : attr.second) {
      std::cout << value.first << ":" << value.second << "; ";
    }
  }
}
