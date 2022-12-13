#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <thread>
#include <vector>

#include "gfd.h"
#include "gfd_validation.h"
#include "graph.h"

namespace algos {
   
void GFDValidation::FitInternal(model::IDatasetStream& data_stream) {}

Eigen::MatrixXi GFDValidation::convert(const Eigen::VectorXi &vec,
                                       const int &cols) {
  Eigen::MatrixXi result = Eigen::MatrixXi::Zero(vec.rows(), cols);
  for (int i = 0; i < vec.rows(); ++i) {
    result(i, vec(i)) = 1;
  }
  return result;
}

bool GFDValidation::isSubgraph(const Pattern &query, const Pattern &graph,
                               const Eigen::MatrixXi &match) {
  Eigen::MatrixXi Q = query.getAdjacencyMatrix();
  Eigen::MatrixXi G = graph.getAdjacencyMatrix();
  if ((match * (match * G).transpose()).transpose().cwiseProduct(Q) != Q) {
    return false;
  }
  Eigen::VectorXi mask = Eigen::VectorXi::Zero(G.rows());
  for (int i = 0; i < G.rows(); ++i) {
    mask(i) = i;
  }
  for (const auto &edge_label : query.getEdges()) {
    std::pair<int, int> edge = edge_label.first;
    Eigen::VectorXi index1 = Eigen::VectorXi::Zero(Q.rows());
    index1(edge.first) = 1;
    int i = mask.cwiseProduct((index1.transpose() * match).transpose()).sum();
    Eigen::VectorXi index2 = Eigen::VectorXi::Zero(Q.rows());
    index2(edge.second) = 1;
    int j = mask.cwiseProduct((index2.transpose() * match).transpose()).sum();
    if (edge_label.second != graph.getEdges().at(std::pair<int, int>(i, j))) {
      return false;
    }
  }
  return true;
}

std::vector<Eigen::VectorXi>
GFDValidation::getCandidateMatches(const Pattern &query, const Pattern &graph,
                                   const std::pair<int, int> &link) {
  std::vector<Eigen::VectorXi> result = {};
  Eigen::MatrixXi Q = query.getAdjacencyMatrix();
  Eigen::MatrixXi G = graph.getAdjacencyMatrix();

  Eigen::MatrixXi graph_degrees = graph.getVertexDegrees();
  Eigen::MatrixXi query_degrees = query.getVertexDegrees();

  std::vector<std::vector<int>> available_vertices = {};
  for (int i = 0; i < Q.rows(); ++i) {
    std::vector<int> temp = {};
    available_vertices.push_back(temp);
  }

  for (int i = 0; i < Q.rows(); ++i) {
    if (i == link.first) {
      available_vertices.at(i).push_back(link.second);
    } else {
      for (int j = 0; j < G.rows(); ++j) {
        if (query.getVertices().at(i) == graph.getVertices().at(j) &&
            query_degrees(i, 0) <= graph_degrees(j, 0) &&
            query_degrees(i, 1) <= graph_degrees(j, 1)) {
          available_vertices.at(i).push_back(j);
        }
      }
    }
  }

  std::vector<Eigen::VectorXi> current = {};
  for (int j = 0; j < available_vertices.at(0).size(); ++j) {
    Eigen::VectorXi temp = Eigen::VectorXi::Zero(Q.rows());
    temp(0) = available_vertices.at(0).at(j);
    current.push_back(temp);
  }
  for (int i = 1; i < available_vertices.size(); ++i) {
    std::vector<Eigen::VectorXi> new_vecs = {};
    for (const Eigen::VectorXi &vec : current) {
      for (int j = 0; j < available_vertices.at(i).size(); ++j) {
        bool contains = false;
        int index = available_vertices.at(i).at(j);
        for (int k = 0; k < i; ++k) {
          if (vec(k) == index) {
            contains = true;
            break;
          }
        }
        if (!contains) {
          Eigen::VectorXi temp = vec;
          temp(i) = index;
          new_vecs.push_back(temp);
        }
      }
    }
    current.clear();
    current = new_vecs;
  }

  for (const Eigen::VectorXi &vec : current) {
    result.push_back(vec);
  }
  return result;
}

std::vector<Eigen::VectorXi>
GFDValidation::getMatches(const Pattern &query, const Pattern &graph,
                          const std::pair<int, int> &link) {
  std::vector<Eigen::VectorXi> result = {};
  std::vector<Eigen::VectorXi> candidates =
      getCandidateMatches(query, graph, link);
  for (const Eigen::VectorXi &match : candidates) {
    if (isSubgraph(query, graph,
                   convert(match, graph.getAdjacencyMatrix().rows()))) {
      result.push_back(match);
    }
  }
  return result;
}

bool GFDValidation::satisfied(const Graph &graph, const Eigen::VectorXi &match,
                              const std::vector<Literal> &literals) {
  for (const Literal &l : literals) {
    auto fst_token = l.first;
    auto snd_token = l.second;
    std::string fst;
    std::string snd;
    if (fst_token.first == -1) {
      fst = fst_token.second;
    } else {
      int index = match(fst_token.first);
      auto attrs = graph.getAttributes().at(index);
      if (attrs.find(fst_token.second) == attrs.end()) {
        return false;
      }
      fst = attrs.at(fst_token.second);
    }
    if (snd_token.first == -1) {
      snd = snd_token.second;
    } else {
      int index = match(snd_token.first);
      auto attrs = graph.getAttributes().at(index);
      if (attrs.find(snd_token.second) == attrs.end()) {
        return false;
      }
      snd = attrs.at(snd_token.second);
    }
    if (fst != snd) {
      return false;
    }
  }
  return true;
}

bool GFDValidation::isSatisfiedLinked(const Graph &graph, const GFD &gfd,
                                      const std::pair<int, int> &link) {
  std::vector<Eigen::VectorXi> matches =
      getMatches(gfd.getPattern(), graph, link);

  for (const Eigen::VectorXi &match : matches) {
    if (!satisfied(graph, match, gfd.getPremises())) {
      continue;
    }

    if (!satisfied(graph, match, gfd.getConclusion())) {
      return false;
    }

    // bool satisfied = true;
    // for (const Literal &l : gfd.getPremises()) {
    //   std::vector<int> vars = l.getVars();
    //   std::vector<std::string> values = l.getValues();
    //   int query_index1 = vars.at(0);
    //   int graph_index1 = match(query_index1);
    //   std::map<std::string, std::string> current_attrs1 =
    //       graph.getAttributes().at(graph_index1);
    //   if (vars.size() == 1) {
    //     // ConstLiteral
    //     if ((current_attrs1.find(values.at(0)) != current_attrs1.end()) &&
    //         (current_attrs1.at(values.at(0)) != values.at(1))) {
    //       // next match
    //       satisfied = false;
    //       break;
    //     }
    //   } else {
    //     // VarLiteral
    //     int query_index2 = vars.at(1);
    //     int graph_index2 = match(query_index2);
    //     std::map<std::string, std::string> current_attrs2 =
    //         graph.getAttributes().at(graph_index2);
    //     if ((current_attrs1.find(values.at(0)) != current_attrs1.end()) &&
    //         (current_attrs2.find(values.at(1)) != current_attrs2.end()) &&
    //         (current_attrs1.at(values.at(0)) !=
    //          current_attrs2.at(values.at(1)))) {
    //       // next match
    //       satisfied = false;
    //       break;
    //     }
    //   }
    // }

    // if (!satisfied) {
    //   break;
    // }

    // for (const Literal &l : gfd.getConclusion()) {
    //   std::vector<int> vars = l.getVars();
    //   std::vector<std::string> values = l.getValues();
    //   int query_index1 = vars.at(0);
    //   int graph_index1 = match(query_index1);
    //   std::map<std::string, std::string> current_attrs1 =
    //       graph.getAttributes().at(graph_index1);
    //   if (vars.size() == 1) {
    //     // ConstLiteral
    //     if ((current_attrs1.find(values.at(0)) == current_attrs1.end()) ||
    //         (current_attrs1.at(values.at(0)) != values.at(1))) {
    //       return false;
    //     }
    //   } else {
    //     // VarLiteral
    //     int query_index2 = vars.at(1);
    //     int graph_index2 = match(query_index2);
    //     std::map<std::string, std::string> current_attrs2 =
    //         graph.getAttributes().at(graph_index2);
    //     if ((current_attrs1.find(values.at(0)) == current_attrs1.end()) ||
    //         (current_attrs2.find(values.at(1)) == current_attrs2.end()) ||
    //         (current_attrs1.at(values.at(0)) !=
    //          current_attrs2.at(values.at(1)))) {
    //       return false;
    //     }
    //   }
    // }
  }
  return true;
}

std::vector<std::vector<int>>
GFDValidation::balanced(const std::vector<int> &weights,
                        const int &processors_num) {
  int m = std::min(processors_num, (int)weights.size());
  std::vector<std::vector<int>> result = {};
  if (weights.begin() == weights.end()) {
    for (int i = 0; i < processors_num; ++i) {
      std::vector<int> temp = {};
      result.push_back(temp);
    }
    return result;
  }
  for (int i = 0; i < m; ++i) {
    // the first value is index
    std::vector<int> temp = {i};
    result.push_back(temp);
  }
  // fill processors initially
  // count optimal
  double optimal = 0;
  int i = 0;
  for (const int &weight : weights) {
    result.at(i++).push_back(weight);
    i = i == m ? 0 : i;
    optimal += weight;
  }
  optimal /= m;
  // sort processors (for convenience)
  for (std::vector<int> &processor : result) {
    std::sort(processor.begin() + 1, processor.end());
  }
  // ALGORITHM
  // 1st step
  std::vector<int> deleted_large = {};
  std::vector<int> deleted_small = {};
  for (std::vector<int> &processor : result) {
    auto border = processor.end();
    for (auto it = --processor.end(); it != processor.begin() + 1; --it) {
      if (*(it - 1) > optimal / 2) {
        deleted_large.push_back(*it);
        border = it;
      } else {
        break;
      }
    }
    processor.erase(border, processor.end());
  }
  // 2nd step
  Eigen::MatrixXi quality = Eigen::MatrixXi::Zero(m, 3);
  for (const std::vector<int> &processor : result) {
    auto last_small = processor.end();
    auto last = processor.end();
    if (*(--processor.end()) > optimal / 2) {
      --last_small;
    }
    if (processor.begin() + 1 == last_small) {
      continue;
    }
    int a = 0;
    int b = 0;
    float sum_small =
        std::accumulate(processor.begin() + 1, last_small, 0, std::plus<int>());
    float sum =
        std::accumulate(processor.begin() + 1, last, 0, std::plus<int>());
    while (sum_small > optimal / 2) {
      ++a;
      --last_small;
      sum_small -= *last_small;
    }
    while (sum > optimal) {
      ++b;
      --last;
      sum -= *last;
    }
    quality(processor.at(0), 0) = a;
    quality(processor.at(0), 1) = b;
    quality(processor.at(0), 2) = a - b;
  }
  // 3rd step
  // sort for convenience
  std::vector<std::vector<int>> small_processors = {};
  std::vector<std::vector<int>> large_processors = {};
  for (const std::vector<int> &processor : result) {
    if (*(--processor.end()) > optimal / 2) {
      large_processors.push_back(processor);
    } else {
      small_processors.push_back(processor);
    }
  }
  auto cGreater = [&quality](std::vector<int> a, std::vector<int> b) {
    return quality(a.at(0), 2) > quality(b.at(0), 2);
  };
  sort(small_processors.begin(), small_processors.end(), cGreater);
  sort(large_processors.begin(), large_processors.end(), cGreater);
  result.clear();
  result.insert(result.end(), small_processors.begin(), small_processors.end());
  result.insert(result.end(), large_processors.begin(), large_processors.end());
  int numOfLarges = large_processors.size() + deleted_large.size();
  // work
  auto border = numOfLarges < m ? result.end() - numOfLarges : result.begin();
  for (auto it = border; it != result.end(); ++it) {
    auto last = it->end();
    if (*(last - 1) > optimal / 2) {
      --last;
    }
    for (auto cur = last - quality(*it->begin(), 0); cur != last; ++cur) {
      deleted_small.push_back(*cur);
    }
    it->erase(last - quality(*it->begin(), 0), last);
  }
  // 4th step
  for (auto it = result.begin(); it != border; ++it) {
    auto last = it->end();
    for (auto cur = last - quality(*it->begin(), 1); cur != last; ++cur) {
      deleted_small.push_back(*cur);
    }
    it->erase(last - quality(*it->begin(), 1), last);
  }
  // 5th step
  i = 0;
  for (const int &weight : deleted_large) {
    if (i < m - large_processors.size()) {
      (result.begin() + i)->push_back(weight);
    } else {
      sort(result.begin(), result.end(),
           [](std::vector<int> a, std::vector<int> b) {
             return std::accumulate(a.begin(), a.end(), 0, std::plus<int>()) <
                    std::accumulate(b.begin(), b.end(), 0, std::plus<int>());
           });
      result.begin()->push_back(weight);
    }
    ++i;
  }
  // 6th step
  for (const int &weight : deleted_small) {
    sort(result.begin(), result.end(),
         [](std::vector<int> a, std::vector<int> b) {
           return std::accumulate(a.begin(), a.end(), 0, std::plus<int>()) <
                  std::accumulate(b.begin(), b.end(), 0, std::plus<int>());
         });
    result.begin()->push_back(weight);
  }
  // delete indices
  for (std::vector<int> &processor : result) {
    processor.erase(processor.begin());
  }
  for (int i = 0; i < processors_num - m; ++i) {
    std::vector<int> empty = {};
    result.push_back(empty);
  }
  return result;
}

Graph GFDValidation::getSubgraph(const Graph &graph, const int &index,
                                 const int &radius, int &out) {
  Eigen::MatrixXi m = graph.getAdjacencyMatrix();
  Eigen::MatrixXi temp = m + m.transpose();
  Eigen::MatrixXi neighbours = Eigen::MatrixXi::Zero(temp.rows(), temp.cols());
  for (int i = 0; i < radius; ++i) {
    Eigen::MatrixXi cur = temp;
    for (int j = 0; j < i; ++j) {
      cur = cur * temp;
    }
    neighbours += cur;
  }
  std::set<int> needed_vertices = {index};

  for (int i = 0; i < neighbours.rows(); ++i) {
    if (neighbours(index, i) != 0) {
      needed_vertices.insert(i);
    }
  }

  std::map<int, std::string> vertices{};
  std::map<int, std::map<std::string, std::string>> attributes{};
  int j = 0;
  std::map<int, int> isomorphism{};
  for (const int &i : needed_vertices) {
    vertices.emplace(j, graph.getVertices().at(i));
    attributes.emplace(j, graph.getAttributes().at(i));
    isomorphism.emplace(i, j++);
  }
  std::map<std::pair<int, int>, std::string> edges{};
  for (const auto &edge_label : graph.getEdges()) {
    std::pair<int, int> edge = edge_label.first;
    if ((needed_vertices.find(edge.first) != needed_vertices.end()) &&
        (needed_vertices.find(edge.second) != needed_vertices.end())) {
      edges.emplace(std::pair<int, int>(isomorphism.at(edge.first),
                                        isomorphism.at(edge.second)),
                    edge_label.second);
    }
  }
  Graph result = Graph(vertices, edges, attributes);
  out = isomorphism.at(index);
  return result;
}

std::vector<int>
GFDValidation::getCandidateVertices(const Pattern &pattern,
                                    const std::string &label,
                                    const std::pair<int, int> &degrees) {
  std::vector<int> result = {};
  Eigen::MatrixXi G = pattern.getAdjacencyMatrix();
  for (const auto &vertex_label : pattern.getVertices()) {
    std::string current_label = vertex_label.second;
    if (current_label != label) {
      continue;
    }
    int index = vertex_label.first;
    if ((G.row(index).sum() >= degrees.first) &&
        (G.col(index).sum() >= degrees.second)) {
      result.push_back(index);
    }
  }
  return result;
}

int GFDValidation::getRadius(const Pattern &pattern, const int &index) {
  Eigen::MatrixXi G = pattern.getAdjacencyMatrix();
  Eigen::MatrixXi m = G + G.transpose();
  Eigen::VectorXi result = Eigen::VectorXi::Zero(G.rows());
  result(index) = 1;
  Eigen::VectorXi mask = Eigen::VectorXi(G.rows());
  mask.setConstant(1);
  mask(index) = 0;
  Eigen::VectorXi prev = Eigen::VectorXi(G.rows());
  prev.setConstant(1);
  int answer = -1;
  while (mask != prev) {
    prev = mask;
    result = (result.transpose() * m).transpose().cwiseProduct(mask);
    // optimize?
    for (int i = 0; i < G.rows(); ++i) {
      if (result(i) != 0) {
        mask(i) = 0;
      }
    }
    ++answer;
  }
  return answer;
}

int GFDValidation::getCenter(const Pattern &pattern, int &out) {
  int min = pattern.getAdjacencyMatrix().rows();
  int result = 0;
  for (const int &index : pattern.getIndices()) {
    int radius = getRadius(pattern, index);
    if (radius <= min) {
      min = radius;
      result = index;
    }
  }
  out = min;
  return result;
}

void GFDValidation::calculateUnsatisfied(
    const std::vector<std::tuple<int, int, int, int>> &messages,
    const std::map<int, GFD> &coded_gfds,
    const std::map<int, Graph> &coded_subgraphs, std::set<int> &out) {
  for (const auto &message : messages) {
    int gfd_index = std::get<0>(message);
    int subgraph_index = std::get<1>(message);
    int center = std::get<2>(message);
    int candidate = std::get<3>(message);
    if (!isSatisfiedLinked(coded_subgraphs.at(subgraph_index),
                           coded_gfds.at(gfd_index),
                           std::pair<int, int>(center, candidate))) {
      out.insert(gfd_index);
    }
  }
}

std::vector<GFD> GFDValidation::getSatisfiedGFDs(const Graph &graph,
                                                 const std::vector<GFD> &gfds,
                                                 const int &m) {
  std::vector<GFD> result = {};
  std::set<int> unsatisfied = {};
  std::map<int, std::vector<std::tuple<int, int, int, int>>> weighted_messages;
  std::vector<int> weights = {};

  std::map<int, GFD> coded_gfds;
  std::map<int, Graph> coded_subgraphs;
  int i = 0;
  int j = 0;
  for (const GFD &gfd : gfds) {
    coded_gfds.emplace(i, gfd);

    int radius = 0;
    int center = getCenter(gfd.getPattern(), radius);
    std::pair<int, int> degrees(gfd.getPattern().getVertexDegrees()(center, 0),
                                gfd.getPattern().getVertexDegrees()(center, 1));
    std::vector<int> candidate_vertices = getCandidateVertices(
        graph, gfd.getPattern().getVertices().at(center), degrees);

    for (const int &candidate : candidate_vertices) {
      int pin = 0;
      Graph subgraph = getSubgraph(graph, candidate, radius, pin);
      coded_subgraphs.emplace(j, subgraph);
      std::tuple<int, int, int, int> temp(i, j, center, pin);
      int weight = subgraph.getSize();
      if (weighted_messages.find(weight) != weighted_messages.end()) {
        weighted_messages.at(weight).push_back(temp);
      } else {
        std::vector<std::tuple<int, int, int, int>> temps = {temp};
        weighted_messages.emplace(weight, temps);
      }
      weights.push_back(weight);
      ++j;
    }
    ++i;
  }

  std::vector<std::vector<int>> balanced_weights = balanced(weights, m);

  std::vector<std::set<int>> answers = {};
  std::vector<std::vector<std::tuple<int, int, int, int>>> groups = {};
  for (int i = 0; i < m; ++i) {
    std::set<int> current = {};
    answers.push_back(current);
    std::vector<std::tuple<int, int, int, int>> messages = {};
    for (int &weight : balanced_weights.at(i)) {
      std::tuple<int, int, int, int> temp =
          *(--weighted_messages.at(weight).end());
      weighted_messages.at(weight).erase(--weighted_messages.at(weight).end());
      messages.push_back(temp);
    }
    groups.push_back(messages);
  }
  std::vector<std::thread> threads = {};
  for (int i = 0; i < m; ++i) {
    std::thread thrd(&GFDValidation::calculateUnsatisfied, this,
                     std::cref(groups.at(i)), std::cref(coded_gfds),
                     std::cref(coded_subgraphs), std::ref(answers.at(i)));
    threads.push_back(std::move(thrd));
  }
  for (std::thread &thrd : threads) {
    if (thrd.joinable()) {
      thrd.join();
    }
  }
  for (const std::set<int> &answer : answers) {
    // optimize?
    for (const int &gfd_index : answer) {
      unsatisfied.insert(gfd_index);
    }
  }
  for (int i = 0; i < gfds.size(); ++i) {
    if (unsatisfied.find(i) == unsatisfied.end()) {
      result.push_back(coded_gfds.at(i));
    }
  }
  return result;
}

unsigned long long GFDValidation::ExecuteInternal() {
  auto start_time = std::chrono::system_clock::now();

  this->result = this->getSatisfiedGFDs(this->graph, this->gfds,
                                        std::thread::hardware_concurrency());

  auto elapsed_milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::system_clock::now() - start_time);
  std::cout << "Satisfied GFDs: " << result.size()
    << "/" << gfds.size() << std::endl;
  return elapsed_milliseconds.count();
}

} // namespace algos
