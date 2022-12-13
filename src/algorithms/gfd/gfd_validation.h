#pragma once
#include <vector>

#include "algorithms/primitive.h"
#include "gfd.h"
#include "graph.h"
#include "parser.h"

namespace algos {

class GFDValidation : public Primitive {
private:
  Graph graph;
  std::vector<GFD> gfds;
  std::vector<GFD> result;

  Eigen::MatrixXi convert(const Eigen::VectorXi &, const int &);
  bool isSubgraph(const Pattern &, const Pattern &, const Eigen::MatrixXi &);
  std::vector<Eigen::VectorXi> getCandidateMatches(const Pattern &,
                                                   const Pattern &,
                                                   const std::pair<int, int> &);
  std::vector<Eigen::VectorXi> getMatches(const Pattern &, const Pattern &,
                                          const std::pair<int, int> &);
  bool satisfied(const Graph &, const Eigen::VectorXi &,
                 const std::vector<Literal> &);
  bool isSatisfiedLinked(const Graph &, const GFD &,
                         const std::pair<int, int> &);
  std::vector<std::vector<int>> balanced(const std::vector<int> &, const int &);
  Graph getSubgraph(const Graph &, const int &, const int &, int &);
  std::vector<int> getCandidateVertices(const Pattern &, const std::string &,
                                        const std::pair<int, int> &);
  int getRadius(const Pattern &, const int &);
  int getCenter(const Pattern &, int &);
  void calculateUnsatisfied(const std::vector<std::tuple<int, int, int, int>> &,
                            const std::map<int, GFD> &,
                            const std::map<int, Graph> &, std::set<int> &);
  std::vector<GFD> getSatisfiedGFDs(const Graph &, const std::vector<GFD> &,
                                    const int &);
  
  void FitInternal(model::IDatasetStream& data_stream);
  unsigned long long ExecuteInternal();
public:
  GFDValidation() : Primitive({}) {};
  GFDValidation(Graph graph_, std::vector<GFD> gfds_)
      : Primitive({}), graph(graph_), gfds(gfds_) { ExecutePrepare(); }
  std::vector<GFD> GFDList() { return this->result; }
};

} // namespace algos
