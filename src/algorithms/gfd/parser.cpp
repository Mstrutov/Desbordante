#include "parser.h"

namespace GFDtools {

std::vector<std::string> Parser::split(std::string str, std::string sep) {
  std::vector<std::string> result = {};
  if (str.empty()) {
    return result;
  }
  size_t pos = 0;
  while ((pos = str.find(sep)) != std::string::npos) {
    result.push_back(str.substr(0, pos));
    str.erase(0, pos + sep.length());
  }
  result.push_back(str);
  return result;
}

std::vector<Literal> Parser::parse_literals(std::istream &stream) {
  std::vector<Literal> result = {};

  std::string line;
  std::getline(stream, line);
  line = line.substr(0, line.find("\r")).substr(0, line.find("\n"));
  auto tokens = Parser::split(line, " ");
  for (auto token : tokens) {
    auto custom_names = Parser::split(token, "=");
    auto names1 = Parser::split(custom_names.at(0), ".");
    int index1 = names1.size() == 1 ? -1 : stoi(names1.at(0));
    std::string name1 = *(--names1.end());
    Token t1(index1, name1);

    auto names2 = Parser::split(custom_names.at(1), ".");
    int index2 = names2.size() == 1 ? -1 : stoi(names2.at(0));
    std::string name2 = *(--names2.end());
    Token t2(index2, name2);

    result.push_back(Literal(t1, t2));
  }

  return result;
}

GFD Parser::parse_gfd(std::istream &stream) {
  std::map<int, std::string> vertices;
  std::map<std::pair<int, int>, std::string> edges;

  std::string line;
  std::getline(stream, line);
  auto sizes = Parser::split(line, " ");
  size_t vertices_num = stoi(sizes.at(0));
  size_t edges_num = stoi(sizes.at(1));
  for (int i = 0; i < vertices_num; ++i) {
    std::getline(stream, line);
    int index = stoi(line.substr(0, line.find("[")));
    line = line.substr(line.find("[") + 1, line.find("]") - line.find("[") - 1);
    auto names = Parser::split(line, "=");
    vertices.emplace(index, names.at(1));
  }
  for (int i = 0; i < edges_num; ++i) {
    std::getline(stream, line);
    auto indices = Parser::split(line.substr(0, line.find("[")), "->");
    int index1 = stoi(indices.at(0));
    int index2 = stoi(indices.at(1));
    line = line.substr(line.find("[") + 1, line.find("]") - line.find("[") - 1);
    auto names = Parser::split(line, "=");
    edges.emplace(std::pair<int, int>(index1, index2), names.at(1));
  }
  std::vector<Literal> premises = Parser::parse_literals(stream);
  std::vector<Literal> conclusion = Parser::parse_literals(stream);

  Pattern pat = Pattern(vertices, edges);
  GFD gfd = GFD(pat, premises, conclusion);
  return gfd;
}

Graph Parser::parse_graph(std::istream &stream) {
  std::map<int, std::string> vertices;
  std::map<std::pair<int, int>, std::string> edges;
  std::map<int, std::map<std::string, std::string>> attributes;

  std::string line;
  std::getline(stream, line);
  auto sizes = Parser::split(line, " ");
  size_t vertices_num = stoi(sizes.at(0));
  size_t edges_num = stoi(sizes.at(1));
  for (int i = 0; i < vertices_num; ++i) {
    std::map<std::string, std::string> attrs;
    std::getline(stream, line);
    int index = stoi(line.substr(0, line.find("[")));
    line = line.substr(line.find("[") + 1, line.find("]") - line.find("[") - 1);
    auto tokens = Parser::split(line, " ");
    for (auto token : tokens) {
      auto names = Parser::split(token, "=");
      std::string fst = names.at(0);
      std::string snd = names.at(1);
      if (fst == "label") {
        vertices.emplace(index, snd);
      } else {
        attrs.emplace(fst, snd);
      }
    }
    attributes.emplace(index, attrs);
  }
  for (int i = 0; i < edges_num; ++i) {
    std::getline(stream, line);
    auto indices = Parser::split(line.substr(0, line.find("[")), "->");
    int index1 = stoi(indices.at(0));
    int index2 = stoi(indices.at(1));
    line = line.substr(line.find("[") + 1, line.find("]") - line.find("[") - 1);
    auto names = Parser::split(line, "=");
    edges.emplace(std::pair<int, int>(index1, index2), names.at(1));
  }
  return Graph(vertices, edges, attributes);
}

}
