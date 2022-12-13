#pragma once
#include <iostream>
#include <map>

#include "gfd.h"
#include "graph.h"

namespace GFDtools {

class Parser {
private:
  static std::vector<std::string> split(std::string, std::string);
  static std::vector<Literal> parse_literals(std::istream &);

public:
  static GFD parse_gfd(std::istream &);
  static Graph parse_graph(std::istream &);
};

}
