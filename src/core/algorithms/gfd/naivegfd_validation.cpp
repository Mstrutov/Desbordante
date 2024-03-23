#include "naivegfd_validation.h"

#include <iostream>
#include <set>

#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "gfd.h"

namespace {

class CheckCallback {
private:
    graph_t const& query_;
    graph_t const& graph_;
    std::vector<Literal> const premises_;
    std::vector<Literal> const conclusion_;
    bool& res_;
    int& amount_;

public:
    CheckCallback(graph_t const& query_, graph_t const& graph_,
                  std::vector<Literal> const& premises_, std::vector<Literal> const& conclusion_,
                  bool& res_, int& amount_)
        : query_(query_),
          graph_(graph_),
          premises_(premises_),
          conclusion_(conclusion_),
          res_(res_),
          amount_(amount_) {}

    template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
    bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) const {
        amount_++;
        auto satisfied = [this, &f](std::vector<Literal> const& literals) {
            for (const Literal& l : literals) {
                auto fst_token = l.first;
                auto snd_token = l.second;
                std::string fst;
                std::string snd;
                if (fst_token.first == -1) {
                    fst = fst_token.second;
                } else {
                    vertex_t v;
                    vertex_t u = boost::vertex(fst_token.first, query_);
                    v = get(f, u);
                    auto attrs = graph_[v].attributes;
                    if (attrs.find(fst_token.second) == attrs.end()) {
                        return false;
                    }
                    fst = attrs.at(fst_token.second);
                }
                if (snd_token.first == -1) {
                    snd = snd_token.second;
                } else {
                    vertex_t v;
                    vertex_t u = boost::vertex(fst_token.first, query_);
                    v = get(f, u);
                    auto attrs = graph_[v].attributes;
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
        };

        if (!satisfied(premises_)) {
            return true;
        }
        if (!satisfied(conclusion_)) {
            res_ = false;
            return false;
        }
        return true;
    }
};

bool Validate(graph_t const& graph, Gfd const& gfd) {
    graph_t pattern = gfd.GetPattern();

    struct VCompare {
        graph_t const& pattern;
        graph_t const& graph;

        bool operator()(vertex_t fr, vertex_t to) const {
            return pattern[fr].attributes.at("label") == graph[to].attributes.at("label");
        }
    } vcompare{pattern, graph};

    struct ECompare {
        graph_t const& pattern;
        graph_t const& graph;

        bool operator()(edge_t fr, edge_t to) const {
            return pattern[fr].label == graph[to].label;
        }
    } ecompare{pattern, graph};

    bool res = true;
    int amount = 0;
    CheckCallback callback(pattern, graph, gfd.GetPremises(), gfd.GetConclusion(), res, amount);

    bool found = boost::vf2_subgraph_iso(
            pattern, graph, callback, get(boost::vertex_index, pattern),
            get(boost::vertex_index, graph), vertex_order_by_mult(pattern), ecompare, vcompare);
    std::cout << "Checked embeddings: " << amount << std::endl;
    if (!found) {
        return true;
    }
    return res;
}

}  // namespace

namespace algos {

std::vector<Gfd> NaiveGfdValidation::GenerateSatisfiedGfds(graph_t const& graph,
                                                           std::vector<Gfd> const& gfds) {
    for (auto& gfd : gfds) {
        if (Validate(graph, gfd)) {
            result_.push_back(gfd);
        }
    }
    return result_;
}

}  // namespace algos
