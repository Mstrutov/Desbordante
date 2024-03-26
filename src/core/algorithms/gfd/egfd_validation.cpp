#include "egfd_validation.h"

#include <iostream>

#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "config/equal_nulls/option.h"
#include "config/names_and_descriptions.h"
#include "config/option_using.h"
#include "config/tabular_data/input_table/option.h"

namespace {

using namespace algos;
using Match = std::vector<std::pair<std::set<vertex_t>::iterator, std::set<vertex_t>::iterator>>;

void FstStepForest(graph_t const& graph, std::map<vertex_t, std::set<vertex_t>>& rooted_subtree,
                   std::map<vertex_t, vertex_t>& children_amount) {
    typename boost::graph_traits<graph_t>::vertex_iterator it, end;
    for (boost::tie(it, end) = vertices(graph); it != end; ++it) {
        if (boost::degree(*it, graph) != 1) {
            continue;
        }
        typename boost::graph_traits<graph_t>::adjacency_iterator adjacency_it =
                boost::adjacent_vertices(*it, graph).first;

        if (rooted_subtree.find(*adjacency_it) != rooted_subtree.end()) {
            rooted_subtree.at(*adjacency_it).insert(*it);
            children_amount[*adjacency_it]++;
        } else {
            std::set<vertex_t> value = {*it};
            rooted_subtree.emplace(*adjacency_it, value);
            children_amount.emplace(*adjacency_it, 1);
        }
    }
}

void BuildForest(graph_t const& graph, std::map<vertex_t, std::set<vertex_t>>& rooted_subtree,
                 std::map<vertex_t, vertex_t>& children_amount) {
    bool changed = true;
    while (changed) {
        changed = false;
        std::map<vertex_t, std::set<vertex_t>> temp;
        std::map<vertex_t, vertex_t> children_temp;
        for (auto const& kv : rooted_subtree) {
            auto desc = kv.first;
            auto children = kv.second;

            if (boost::degree(desc, graph) == (children_amount.at(desc) + 1)) {
                changed = true;
                typename boost::graph_traits<graph_t>::adjacency_iterator adjacency_it,
                        adjacency_end;
                boost::tie(adjacency_it, adjacency_end) = boost::adjacent_vertices(desc, graph);
                for (; adjacency_it != adjacency_end; ++adjacency_it) {
                    if (children.find(*adjacency_it) != children.end()) {
                        continue;
                    }
                    if (temp.find(*adjacency_it) != temp.end()) {
                        std::set<vertex_t> value = {};
                        auto current = temp.at(*adjacency_it);
                        std::set_union(children.begin(), children.end(), current.begin(),
                                       current.end(), std::inserter(value, value.begin()));
                        value.insert(desc);
                        temp[*adjacency_it] = value;
                        children_temp[*adjacency_it]++;
                    } else {
                        std::set<vertex_t> value = children;
                        value.insert(desc);
                        temp.emplace(*adjacency_it, value);
                        children_temp.emplace(*adjacency_it, 1);
                    }
                    break;
                }
            } else {
                temp.emplace(desc, children);
                children_temp.emplace(desc, children_amount.at(desc));
            }
        }
        rooted_subtree = temp;
        children_amount = children_temp;
    }
}

void CfDecompose(graph_t const& graph, std::set<vertex_t>& core,
                 std::vector<std::set<vertex_t>>& forest) {
    if (boost::num_vertices(graph) == (boost::num_edges(graph) + 1)) {
        typename boost::graph_traits<graph_t>::vertex_iterator it, end;
        for (boost::tie(it, end) = vertices(graph); it != end; ++it) {
            core.insert(*it);
        }
        return;
    }

    std::map<vertex_t, std::set<vertex_t>> rooted_subtree;
    std::map<vertex_t, vertex_t> children_amount;
    FstStepForest(graph, rooted_subtree, children_amount);
    BuildForest(graph, rooted_subtree, children_amount);

    std::set<vertex_t> not_core_indices = {};
    for (auto const& kv : rooted_subtree) {
        for (int child : kv.second) {
            not_core_indices.insert(child);
        }
    }
    typename boost::graph_traits<graph_t>::vertex_iterator it, end;
    for (boost::tie(it, end) = vertices(graph); it != end; ++it) {
        if (not_core_indices.find(*it) == not_core_indices.end()) {
            core.insert(*it);
        }
    }

    for (auto const& kv : rooted_subtree) {
        std::set<vertex_t> indices(kv.second);
        indices.insert(kv.first);
        forest.push_back(indices);
    }
}

int Mnd(graph_t const& graph, vertex_t const& v) {
    typename boost::graph_traits<graph_t>::adjacency_iterator adjacency_it, adjacency_end;
    boost::tie(adjacency_it, adjacency_end) = boost::adjacent_vertices(v, graph);
    std::size_t result = 0;
    for (; adjacency_it != adjacency_end; ++adjacency_it) {
        if (result < boost::degree(*adjacency_it, graph)) {
            result = boost::degree(*adjacency_it, graph);
        }
    }
    return result;
}

void CountLabelDegrees(graph_t const& graph, vertex_t const& v,
                       std::map<std::string, vertex_t>& result) {
    typename boost::graph_traits<graph_t>::adjacency_iterator adjacency_it, adjacency_end;
    boost::tie(adjacency_it, adjacency_end) = boost::adjacent_vertices(v, graph);
    for (; adjacency_it != adjacency_end; ++adjacency_it) {
        if (result.find(graph[*adjacency_it].attributes.at("label")) != result.end()) {
            result[graph[*adjacency_it].attributes.at("label")]++;
        } else {
            result.emplace(graph[*adjacency_it].attributes.at("label"), 1);
        }
    }
}

bool CandVerify(graph_t const& graph, vertex_t const& v, graph_t const& query, vertex_t const& u) {
    if (Mnd(graph, v) < Mnd(query, u)) {
        return false;
    }
    std::map<std::string, vertex_t> graph_label_degrees;
    CountLabelDegrees(graph, v, graph_label_degrees);
    std::map<std::string, vertex_t> query_label_degrees;
    CountLabelDegrees(query, u, query_label_degrees);

    for (auto const& label_degree : query_label_degrees) {
        std::string const& label = label_degree.first;
        std::size_t const& degree = label_degree.second;
        if (graph_label_degrees.find(label) == graph_label_degrees.end() ||
            graph_label_degrees.at(label) < degree) {
            return false;
        }
    }
    return true;
}

void SortComplexity(std::vector<vertex_t>& order, graph_t const& graph, graph_t const& query,
                    std::map<std::string, std::set<vertex_t>> const& label_classes) {
    auto cmp_complexity = [&graph, &query, &label_classes](vertex_t const& a, vertex_t const& b) {
        std::size_t a_degree = boost::degree(a, query);
        int an = 0;
        for (const vertex_t& e : label_classes.at(query[a].attributes.at("label"))) {
            if (boost::degree(e, graph) >= a_degree) {
                an++;
            }
        }

        std::size_t b_degree = boost::degree(b, query);
        int bn = 0;
        for (const vertex_t& e : label_classes.at(query[b].attributes.at("label"))) {
            if (boost::degree(e, graph) >= b_degree) {
                bn++;
            }
        }
        return an / a_degree < bn / b_degree;
    };
    std::sort(order.begin(), order.end(), cmp_complexity);
}

void SortAccurateComplexity(std::vector<vertex_t>& order, graph_t const& graph,
                            graph_t const& query,
                            std::map<std::string, std::set<vertex_t>> const& label_classes) {
    int top = std::min(int(order.size()), 3);
    auto cmp_accurate_complexity = [&graph, &query, &label_classes](vertex_t const& a,
                                                                    vertex_t const& b) {
        int a_degree = boost::degree(a, query);
        int an = 0;
        for (const vertex_t& e : label_classes.at(query[a].attributes.at("label"))) {
            if (CandVerify(graph, e, query, a)) {
                an++;
            }
        }

        int b_degree = boost::degree(b, query);
        int bn = 0;
        for (const vertex_t& e : label_classes.at(query[b].attributes.at("label"))) {
            if (CandVerify(graph, e, query, b)) {
                bn++;
            }
        }
        return an / a_degree < bn / b_degree;
    };
    std::sort(order.begin(), std::next(order.begin(), top), cmp_accurate_complexity);
}

int GetRoot(graph_t const& graph, graph_t const& query, std::set<vertex_t> const& core) {
    std::map<std::string, std::set<vertex_t>> label_classes;
    typename boost::graph_traits<graph_t>::vertex_iterator it, end;
    for (boost::tie(it, end) = vertices(graph); it != end; ++it) {
        if (label_classes.find(graph[*it].attributes.at("label")) != label_classes.end()) {
            label_classes[graph[*it].attributes.at("label")].insert(*it);
        } else {
            std::set<vertex_t> value = {*it};
            label_classes.emplace(graph[*it].attributes.at("label"), value);
        }
    }
    std::vector<vertex_t> order(core.begin(), core.end());

    SortComplexity(order, graph, query, label_classes);
    SortAccurateComplexity(order, graph, query, label_classes);
    return *order.begin();
}

void MakeLevels(graph_t const& query, vertex_t const& root, std::vector<std::set<vertex_t>>& levels,
                std::map<vertex_t, vertex_t>& parent) {
    std::set<vertex_t> current = {root};
    std::set<vertex_t> marked = {root};
    while (!current.empty()) {
        levels.push_back(current);
        std::set<vertex_t> next = {};
        for (vertex_t const& vertex : current) {
            typename boost::graph_traits<graph_t>::adjacency_iterator adjacency_it, adjacency_end;
            boost::tie(adjacency_it, adjacency_end) = boost::adjacent_vertices(vertex, query);
            for (; adjacency_it != adjacency_end; ++adjacency_it) {
                if (marked.find(*adjacency_it) == marked.end()) {
                    marked.insert(*adjacency_it);
                    next.insert(*adjacency_it);
                    parent.emplace(*adjacency_it, vertex);
                }
            }
        }
        current = next;
    }
}

void MakeNte(graph_t const& query, std::vector<std::set<vertex_t>>& levels,
             std::map<vertex_t, vertex_t>& parent, std::set<edge_t>& nte, std::set<edge_t>& snte) {
    typename boost::graph_traits<graph_t>::edge_iterator it_edge, end_edge;
    for (boost::tie(it_edge, end_edge) = edges(query); it_edge != end_edge; ++it_edge) {
        vertex_t origin = boost::source(*it_edge, query);
        vertex_t finish = boost::target(*it_edge, query);
        if ((parent.find(origin) != parent.end()) && (parent.find(finish) != parent.end()) &&
            (parent.at(origin) != finish) && (parent.at(finish) != origin)) {
            int origin_level = 0;
            int finish_level = 0;
            for (std::size_t i = 0; i < levels.size(); ++i) {
                if (levels.at(i).find(origin) != levels.at(i).end()) {
                    origin_level = i;
                }
                if (levels.at(i).find(finish) != levels.at(i).end()) {
                    finish_level = i;
                }
            }
            if (origin_level == finish_level) {
                snte.insert(*it_edge);
            }
            nte.insert(*it_edge);
        }
    }
}

void BfsTree(graph_t const& query, vertex_t const& root, std::vector<std::set<vertex_t>>& levels,
             std::map<vertex_t, vertex_t>& parent, std::set<edge_t>& nte, std::set<edge_t>& snte) {
    MakeLevels(query, root, levels, parent);
    MakeNte(query, levels, parent, nte, snte);
}

void DirectConstruction(std::set<vertex_t> const& lev, graph_t const& graph, graph_t const& query,
                        std::map<vertex_t, std::set<vertex_t>>& candidates,
                        std::map<vertex_t, int>& cnts,
                        std::map<vertex_t, std::set<vertex_t>>& unvisited_neighbours,
                        std::set<edge_t> const& snte, std::set<vertex_t>& visited) {
    for (vertex_t const& u : lev) {
        int cnt = 0;
        typename boost::graph_traits<graph_t>::adjacency_iterator adjacency_it, adjacency_end;
        boost::tie(adjacency_it, adjacency_end) = boost::adjacent_vertices(u, query);
        for (; adjacency_it != adjacency_end; ++adjacency_it) {
            if (visited.find(query[*adjacency_it].node_id) == visited.end() &&
                snte.find(boost::edge(*adjacency_it, u, query).first) != snte.end()) {
                if (unvisited_neighbours.find(u) != unvisited_neighbours.end()) {
                    unvisited_neighbours.at(u).insert(*adjacency_it);
                } else {
                    std::set<vertex_t> value = {*adjacency_it};
                    unvisited_neighbours.emplace(u, value);
                }
            } else if (visited.find(*adjacency_it) != visited.end()) {
                for (vertex_t const& v : candidates.at(*adjacency_it)) {
                    typename boost::graph_traits<graph_t>::adjacency_iterator g_adj_it, g_adj_end;
                    boost::tie(g_adj_it, g_adj_end) = boost::adjacent_vertices(v, graph);
                    for (; g_adj_it != g_adj_end; ++g_adj_it) {
                        if (graph[*g_adj_it].attributes.at("label") ==
                                    query[u].attributes.at("label") &&
                            boost::degree(*g_adj_it, graph) >= boost::degree(u, query)) {
                            if (cnts.find(*g_adj_it) == cnts.end()) {
                                if (cnt == 0) {
                                    cnts.emplace(*g_adj_it, 1);
                                }
                            } else {
                                if (cnts.at(*g_adj_it) == cnt) {
                                    cnts[*g_adj_it]++;
                                }
                            }
                        }
                    }
                }
                cnt++;
            }
        }
        typename boost::graph_traits<graph_t>::vertex_iterator g_it, g_end;
        for (boost::tie(g_it, g_end) = vertices(graph); g_it != g_end; ++g_it) {
            if (((cnts.find(*g_it) == cnts.end()) && (cnt == 0)) ||
                ((cnts.find(*g_it) != cnts.end()) && (cnts.at(*g_it) == cnt))) {
                if (CandVerify(graph, *g_it, query, u)) {
                    candidates.at(u).insert(*g_it);
                }
            }
        }
        visited.insert(u);
        cnts.clear();
    }
}

void ReverseConstruction(std::set<vertex_t> const& lev, graph_t const& graph, graph_t const& query,
                         std::map<vertex_t, std::set<vertex_t>>& candidates,
                         std::map<vertex_t, int>& cnts,
                         std::map<vertex_t, std::set<vertex_t>>& unvisited_neighbours) {
    for (std::set<vertex_t>::iterator j = --lev.end(); j != std::next(lev.begin(), -1); --j) {
        vertex_t u = *j;
        int cnt = 0;
        if (unvisited_neighbours.find(u) != unvisited_neighbours.end()) {
            for (vertex_t const& un : unvisited_neighbours.at(u)) {
                for (vertex_t const& v : candidates.at(un)) {
                    typename boost::graph_traits<graph_t>::adjacency_iterator g_adj_it, g_adj_end;
                    boost::tie(g_adj_it, g_adj_end) =
                            boost::adjacent_vertices(boost::vertex(v, graph), graph);
                    for (; g_adj_it != g_adj_end; ++g_adj_it) {
                        if (graph[*g_adj_it].attributes.at("label") ==
                                    query[u].attributes.at("label") &&
                            boost::degree(*g_adj_it, graph) >= boost::degree(u, query)) {
                            if (cnts.find(*g_adj_it) == cnts.end()) {
                                if (cnt == 0) {
                                    cnts.emplace(*g_adj_it, 1);
                                }
                            } else {
                                if (cnts.at(*g_adj_it) == cnt) {
                                    cnts[*g_adj_it]++;
                                }
                            }
                        }
                    }
                }
                cnt++;
            }
        }

        std::set<vertex_t> to_delete = {};
        for (vertex_t const& v : candidates.at(u)) {
            if (!(((cnts.find(v) == cnts.end()) && (cnt == 0)) ||
                  ((cnts.find(v) != cnts.end()) && (cnts.at(v) == cnt)))) {
                to_delete.insert(v);
            }
        }
        for (vertex_t const& d : to_delete) {
            candidates.at(u).erase(d);
        }
        cnts.clear();
    }
}

void FinalConstruction(std::set<vertex_t> const& lev, CPI& cpi, graph_t const& graph,
                       graph_t const& query, std::map<vertex_t, vertex_t> const& parent,
                       std::map<vertex_t, std::set<vertex_t>>& candidates) {
    for (vertex_t const& u : lev) {
        vertex_t up = parent.at(u);
        for (vertex_t const& vp : candidates.at(up)) {
            typename boost::graph_traits<graph_t>::adjacency_iterator g_adj_it, g_adj_end;
            boost::tie(g_adj_it, g_adj_end) = boost::adjacent_vertices(vp, graph);
            for (; g_adj_it != g_adj_end; ++g_adj_it) {
                if (graph[*g_adj_it].attributes.at("label") == query[u].attributes.at("label") &&
                    boost::degree(*g_adj_it, graph) >= boost::degree(u, query) &&
                    candidates.at(u).find(*g_adj_it) != candidates.at(u).end() &&
                    query[boost::edge(up, u, query).first].label ==
                            graph[boost::edge(vp, *g_adj_it, graph).first].label) {
                    std::pair<vertex_t, vertex_t> cpi_edge(up, u);
                    if (cpi.find(cpi_edge) != cpi.end()) {
                        if (cpi.at(cpi_edge).find(vp) != cpi.at(cpi_edge).end()) {
                            cpi.at(cpi_edge).at(vp).insert(*g_adj_it);
                        } else {
                            std::set<vertex_t> value = {*g_adj_it};
                            cpi.at(cpi_edge).emplace(vp, value);
                        }
                    } else {
                        std::map<vertex_t, std::set<vertex_t>> edge_map;
                        std::set<vertex_t> value = {*g_adj_it};
                        edge_map.emplace(vp, value);
                        cpi.emplace(cpi_edge, edge_map);
                    }
                }
            }
        }
    }
}

void TopDownConstruct(CPI& cpi, graph_t const& graph, graph_t const& query,
                      std::vector<std::set<vertex_t>> const& levels,
                      std::map<vertex_t, vertex_t> const& parent,
                      std::map<vertex_t, std::set<vertex_t>>& candidates,
                      std::set<edge_t> const& snte) {
    vertex_t root = *levels.at(0).begin();
    typename boost::graph_traits<graph_t>::vertex_iterator it, end;
    for (boost::tie(it, end) = vertices(query); it != end; ++it) {
        std::set<vertex_t> empty = {};
        candidates.emplace(*it, empty);
    }

    for (boost::tie(it, end) = vertices(graph); it != end; ++it) {
        if (graph[*it].attributes.at("label") == query[root].attributes.at("label") &&
            boost::degree(*it, graph) >= boost::degree(root, query) &&
            CandVerify(graph, *it, query, root)) {
            candidates.at(root).insert(*it);
        }
    }
    std::set<vertex_t> visited = {root};
    std::map<vertex_t, std::set<vertex_t>> unvisited_neighbours;
    std::map<vertex_t, int> cnts;

    std::vector<std::set<vertex_t>>::const_iterator i = std::next(levels.cbegin());
    for (; i != levels.cend(); ++i) {
        std::set<vertex_t> lev = *i;
        DirectConstruction(lev, graph, query, candidates, cnts, unvisited_neighbours, snte,
                           visited);
        ReverseConstruction(lev, graph, query, candidates, cnts, unvisited_neighbours);
        FinalConstruction(lev, cpi, graph, query, parent, candidates);
    }
}

void InitialRefinement(vertex_t const& u, graph_t const& graph, graph_t const& query,
                       std::map<vertex_t, vertex_t> const& parent,
                       std::map<vertex_t, std::set<vertex_t>>& candidates,
                       std::map<vertex_t, int>& cnts, int& cnt) {
    typename boost::graph_traits<graph_t>::adjacency_iterator q_adj_it, q_adj_end;
    boost::tie(q_adj_it, q_adj_end) = boost::adjacent_vertices(u, query);
    for (; q_adj_it != q_adj_end; ++q_adj_it) {
        if ((parent.find(*q_adj_it) != parent.end()) && (parent.at(*q_adj_it) == u)) {
            for (vertex_t const& v : candidates.at(*q_adj_it)) {
                typename boost::graph_traits<graph_t>::adjacency_iterator g_adj_it, g_adj_end;
                boost::tie(g_adj_it, g_adj_end) = boost::adjacent_vertices(v, graph);
                for (; g_adj_it != g_adj_end; ++g_adj_it) {
                    if (graph[*g_adj_it].attributes.at("label") ==
                                query[u].attributes.at("label") &&
                        boost::degree(*g_adj_it, graph) >= boost::degree(u, query)) {
                        if (cnts.find(*g_adj_it) == cnts.end()) {
                            if (cnt == 0) {
                                cnts.emplace(*g_adj_it, 1);
                            }
                        } else {
                            if (cnts.at(*g_adj_it) == cnt) {
                                cnts[*g_adj_it]++;
                            }
                        }
                    }
                }
            }
            cnt++;
        }
    }
}

void OddDeletion(vertex_t const& u, CPI& cpi, std::map<vertex_t, std::set<vertex_t>>& candidates,
                 std::map<vertex_t, int>& cnts, int& cnt) {
    std::set<vertex_t> to_delete = {};
    for (vertex_t const& v : candidates.at(u)) {
        if (!(((cnts.find(v) == cnts.end()) && (cnt == 0)) ||
              ((cnts.find(v) != cnts.end()) && (cnts.at(v) == cnt)))) {
            to_delete.insert(v);
        }
    }
    for (vertex_t const& d : to_delete) {
        candidates.at(u).erase(d);
        for (auto const& e : cpi) {
            if (e.second.find(d) != e.second.end()) {
                cpi.at(e.first).erase(d);
            }
        }
    }
    cnts.clear();
}

void FinalRefinement(vertex_t const& u, CPI& cpi, graph_t const& query,
                     std::map<vertex_t, vertex_t> const& parent,
                     std::map<vertex_t, std::set<vertex_t>>& candidates) {
    for (vertex_t const& v : candidates.at(u)) {
        typename boost::graph_traits<graph_t>::adjacency_iterator q_adj_it, q_adj_end;
        boost::tie(q_adj_it, q_adj_end) = boost::adjacent_vertices(u, query);
        for (; q_adj_it != q_adj_end; ++q_adj_it) {
            vertex_t u2 = *q_adj_it;
            if ((parent.find(u2) != parent.end()) && (parent.at(u2) == u)) {
                std::pair<vertex_t, vertex_t> cpi_edge(u, u2);
                for (vertex_t const& v2 : cpi.at(cpi_edge).at(v)) {
                    if (candidates.at(u2).find(v2) == candidates.at(u2).end()) {
                        cpi.at(cpi_edge).at(v).erase(v2);
                    }
                }
            }
        }
    }
}

void BottomUpRefinement(CPI& cpi, graph_t const& graph, graph_t const& query,
                        std::vector<std::set<vertex_t>> const& levels,
                        std::map<vertex_t, vertex_t> const& parent,
                        std::map<vertex_t, std::set<vertex_t>>& candidates) {
    std::map<vertex_t, int> cnts;

    std::vector<std::set<vertex_t>>::const_iterator lev_it;
    for (lev_it = --levels.cend(); lev_it != std::next(levels.begin(), -1); --lev_it) {
        for (vertex_t const& u : *lev_it) {
            int cnt = 0;
            InitialRefinement(u, graph, query, parent, candidates, cnts, cnt);
            OddDeletion(u, cpi, candidates, cnts, cnt);
            FinalRefinement(u, cpi, query, parent, candidates);
        }
    }
}

int NumOfEmbeddings(const CPI& cpi, std::vector<vertex_t> const& path, vertex_t const& origin) {
    std::map<std::pair<vertex_t, vertex_t>, int> result;
    std::pair<vertex_t, vertex_t> edge(*(std::next(path.end(), -2)), *(std::next(path.end(), -1)));
    for (auto& vert_cans : cpi.at(edge)) {
        for (vertex_t const& can : vert_cans.second) {
            std::pair<vertex_t, vertex_t> key(*(std::next(path.end(), -1)), can);
            result.emplace(key, 1);
        }
    }

    for (int i = path.size() - 1; path.at(i) != origin; --i) {
        std::map<std::pair<vertex_t, vertex_t>, int> new_result;
        std::pair<vertex_t, vertex_t> cur_edge(path.at(i - 1), path.at(i));
        for (auto& vert_cans : cpi.at(cur_edge)) {
            for (vertex_t const& can : vert_cans.second) {
                std::pair<vertex_t, vertex_t> key(path.at(i - 1), vert_cans.first);
                std::pair<vertex_t, vertex_t> counted(path.at(i), can);
                if (new_result.find(key) != new_result.end()) {
                    new_result[key] += result.at(counted);
                } else {
                    new_result.emplace(key, result.at(counted));
                }
            }
        }
        result = new_result;
    }
    int answer = 0;
    for (auto& values : result) {
        answer += values.second;
    }
    return answer;
}

int CandidatesCardinality(const CPI& cpi, vertex_t const& u) {
    for (auto& edge_cans : cpi) {
        auto edge = edge_cans.first;
        if (edge.first == u) {
            return edge_cans.second.size();
        }
    }
    return 1;
}

void BuildOptimalSeq(const CPI& cpi, std::vector<std::vector<vertex_t>> const& paths_origin,
                     std::vector<vertex_t> const& NTs,
                     std::vector<std::vector<vertex_t>> const& paths, std::vector<vertex_t>& pi) {
    auto cmp = [&cpi, &paths_origin, &NTs](std::vector<vertex_t> const& a,
                                           std::vector<vertex_t> const& b) {
        int nta = NTs.at(std::find(paths_origin.begin(), paths_origin.end(), a) -
                         paths_origin.begin());
        int ntb = NTs.at(std::find(paths_origin.begin(), paths_origin.end(), b) -
                         paths_origin.begin());
        if (nta == 0) {
            return false;
        }
        if (ntb == 0) {
            return true;
        }
        return NumOfEmbeddings(cpi, a, *a.begin()) / nta <
               NumOfEmbeddings(cpi, b, *b.begin()) / ntb;
    };
    pi = *std::min_element(paths.begin(), paths.end(), cmp);
}

void BuildAccurateOptimalSeq(const CPI& cpi, std::vector<std::vector<vertex_t>> const& paths_origin,
                             std::vector<vertex_t> const& origins,
                             std::vector<std::vector<vertex_t>> const& paths,
                             std::vector<vertex_t>& pi) {
    auto cmp = [&cpi, &paths_origin, &origins](std::vector<vertex_t> const& a,
                                               std::vector<vertex_t> const& b) {
        int a_origin = origins.at(std::find(paths_origin.begin(), paths_origin.end(), a) -
                                  paths_origin.begin());
        int b_origin = origins.at(std::find(paths_origin.begin(), paths_origin.end(), b) -
                                  paths_origin.begin());
        return NumOfEmbeddings(cpi, a, a_origin) / CandidatesCardinality(cpi, a_origin) <
               NumOfEmbeddings(cpi, b, b_origin) / CandidatesCardinality(cpi, b_origin);
    };
    pi = *std::min_element(paths.begin(), paths.end(), cmp);
}

std::vector<vertex_t> MatchingOrder(const CPI& cpi,
                                    std::vector<std::vector<vertex_t>> const& paths_origin,
                                    std::vector<vertex_t> const& NTs) {
    std::vector<std::vector<vertex_t>> paths;
    std::copy(paths_origin.begin(), paths_origin.end(), std::back_inserter(paths));

    std::vector<vertex_t> pi;
    BuildOptimalSeq(cpi, paths_origin, NTs, paths, pi);
    std::vector<vertex_t> origins;
    for (auto& path : paths_origin) {
        if (path == pi) {
            origins.push_back(*path.begin());
            continue;
        }
        std::vector<vertex_t>::iterator pi_it = pi.begin();
        std::vector<vertex_t>::const_iterator path_it = path.cbegin();
        for (; *pi_it == *path_it; ++pi_it, ++path_it) {
        }
        origins.push_back(*(--pi_it));
    }
    paths.erase(std::remove(paths.begin(), paths.end(), pi), paths.end());

    std::vector<vertex_t> seq;
    std::copy(pi.begin(), pi.end(), std::back_inserter(seq));
    while (!paths.empty()) {
        std::vector<vertex_t> pi_new;
        BuildAccurateOptimalSeq(cpi, paths_origin, origins, paths, pi_new);
        std::vector<vertex_t>::iterator pi_it = pi_new.begin();
        std::vector<vertex_t>::iterator seq_it = seq.begin();
        for (; *pi_it == *seq_it; ++pi_it, ++seq_it) {
        }
        seq.insert(seq_it, pi_it, pi_new.end());
        paths.erase(std::remove(paths.begin(), paths.end(), pi_new), paths.end());
    }
    return seq;
}

bool ValidateNt(graph_t const& graph, vertex_t const& v, graph_t const& query, vertex_t const& u,
                std::vector<vertex_t> const& seq, std::map<vertex_t, vertex_t> const& parent,
                Match match) {
    int index = std::find(seq.begin(), seq.end(), u) - seq.begin();
    for (int i = 0; i < index; ++i) {
        if ((seq.at(i) != parent.at(u)) && boost::edge(seq.at(i), u, query).second) {
            if (!boost::edge(*match.at(i).first, v, graph).second ||
                graph[boost::edge(*match.at(i).first, v, graph).first].label !=
                        query[boost::edge(seq.at(i), u, query).first].label) {
                return false;
            }
        }
    }
    return true;
}

std::vector<std::vector<vertex_t>> GetPaths(std::set<vertex_t> const& indices,
                                            std::map<vertex_t, vertex_t> const& parent_) {
    std::vector<std::vector<vertex_t>> result = {};
    std::map<vertex_t, vertex_t> parent(parent_);
    std::set<vertex_t> to_delete = {};
    for (auto& link : parent) {
        if (indices.find(link.first) == indices.end()) {
            to_delete.insert(link.first);
        }
    }
    for (auto& index : to_delete) {
        parent.erase(index);
    }
    std::set<vertex_t> keys = {};
    std::set<vertex_t> values = {};
    for (auto& kv : parent) {
        keys.insert(kv.first);
        values.insert(kv.second);
    }
    std::set<vertex_t> leaves = {};
    std::set_difference(keys.begin(), keys.end(), values.begin(), values.end(),
                        std::inserter(leaves, leaves.begin()));
    for (vertex_t const& leaf : leaves) {
        std::vector<vertex_t> path = {leaf};
        vertex_t cur = leaf;
        while ((parent.find(cur) != parent.end()) &&
               (indices.find(parent.at(cur)) != indices.end())) {
            cur = parent.at(cur);
            path.push_back(cur);
        }
        std::reverse(path.begin(), path.end());
        result.push_back(path);
    }
    return result;
}

bool Visited(Match const& match, vertex_t const& v, std::size_t const& index) {
    for (std::size_t i = 0; i < match.size(); ++i) {
        if (i != index) {
            if ((match.at(i).first != match.at(i).second) && (*match.at(i).first == v)) {
                return true;
            }
        }
    }
    return false;
}

bool Satisfied(graph_t const& graph, graph_t const& query, std::vector<vertex_t> const& seq,
               Match const& match, std::vector<Literal> const& literals) {
    for (Literal const& l : literals) {
        auto fst_token = l.first;
        auto snd_token = l.second;
        std::string fst;
        std::string snd;
        if (fst_token.first == -1) {
            fst = fst_token.second;
        } else {
            vertex_t v;
            vertex_t u = boost::vertex(fst_token.first, query);
            int index = std::find(seq.begin(), seq.end(), u) - seq.begin();
            v = *match.at(index).first;
            auto attrs = graph[v].attributes;
            if (attrs.find(fst_token.second) == attrs.end()) {
                return false;
            }
            fst = attrs.at(fst_token.second);
        }
        if (snd_token.first == -1) {
            snd = snd_token.second;
        } else {
            vertex_t v;
            vertex_t u = boost::vertex(snd_token.first, query);
            int index = std::find(seq.begin(), seq.end(), u) - seq.begin();
            v = *match.at(index).first;
            auto attrs = graph[v].attributes;
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

void FullNTs(std::vector<std::vector<vertex_t>> const& paths, std::set<edge_t> const& nte,
             graph_t const& query, std::vector<vertex_t>& NTs) {
    for (auto& path : paths) {
        int nt = 0;
        for (auto& desc : nte) {
            vertex_t source = boost::source(desc, query);
            vertex_t target = boost::target(desc, query);
            if ((std::find(path.begin(), path.end(), source) != path.end()) ||
                (std::find(path.begin(), path.end(), target) != path.end())) {
                nt++;
            }
        }
        NTs.push_back(nt);
    }
}

void CompleteSeq(CPI& cpi, std::vector<std::set<vertex_t>> const& forest,
                 std::map<vertex_t, vertex_t> const& parent, graph_t const& query,
                 std::set<edge_t> const& nte, std::vector<vertex_t>& seq) {
    for (auto& tree : forest) {
        std::vector<std::vector<vertex_t>> tree_paths = GetPaths(tree, parent);

        std::vector<vertex_t> tree_nts = {};
        for (auto& path : tree_paths) {
            int nt = 0;
            for (auto& desc : nte) {
                vertex_t source = boost::source(desc, query);
                vertex_t target = boost::target(desc, query);
                if ((std::find(path.begin(), path.end(), source) != path.end()) ||
                    (std::find(path.begin(), path.end(), target) != path.end())) {
                    nt++;
                }
            }
            tree_nts.push_back(nt);
        }
        std::vector<vertex_t> tree_seq = MatchingOrder(cpi, tree_paths, tree_nts);

        seq.insert(seq.end(), ++tree_seq.begin(), tree_seq.end());
    }
}

bool FullMatch(CPI& cpi, Match& match, std::set<vertex_t> const& root_candidates,
               std::set<vertex_t> const& core, std::vector<vertex_t> const& seq,
               std::map<vertex_t, vertex_t> const& parent, graph_t const& graph,
               graph_t const& query) {
    match.push_back({root_candidates.begin(), root_candidates.end()});
    for (std::size_t i = 1; i < core.size(); ++i) {
        std::pair<vertex_t, vertex_t> edge(parent.at(seq.at(i)), seq.at(i));
        int index = std::find(seq.begin(), seq.end(), parent.at(seq.at(i))) - seq.begin();
        std::pair<std::set<vertex_t>::iterator, std::set<vertex_t>::iterator> its(
                cpi.at(edge).at(*match.at(index).first).begin(),
                cpi.at(edge).at(*match.at(index).first).end());
        match.push_back(its);

        while ((match.at(i).first != match.at(i).second) &&
               (Visited(match, *match.at(i).first, i) ||
                !ValidateNt(graph, *match.at(i).first, query, seq.at(i), seq, parent, match))) {
            match.at(i).first++;
        }
        if (match.at(i).first == match.at(i).second) {
            // std::cout << "Trivially satisfied" << std::endl;
            return true;
        }
    }
    for (std::size_t i = core.size(); i < seq.size(); ++i) {
        match.push_back({root_candidates.end(), root_candidates.end()});
    }
    return false;
}

void IncrementMatch(int& i, const CPI& cpi, Match& match,
                    std::map<vertex_t, vertex_t> const& parent, std::set<vertex_t> const& core,
                    std::vector<vertex_t> const& seq, graph_t const& graph, graph_t const& query) {
    while ((i != static_cast<int>(core.size())) && (i != -1)) {
        if (match.at(i).first == match.at(i).second) {
            std::pair<vertex_t, vertex_t> edge(parent.at(seq.at(i)), seq.at(i));
            std::size_t index =
                    std::find(seq.begin(), seq.end(), parent.at(seq.at(i))) - seq.begin();
            std::pair<std::set<vertex_t>::iterator, std::set<vertex_t>::iterator> its(
                    cpi.at(edge).at(*match.at(index).first).begin(),
                    cpi.at(edge).at(*match.at(index).first).end());
            match[i] = its;
        } else {
            match.at(i).first++;
        }

        while ((match.at(i).first != match.at(i).second) &&
               (Visited(match, *match.at(i).first, i) ||
                !ValidateNt(graph, *match.at(i).first, query, seq.at(i), seq, parent, match))) {
            match.at(i).first++;
        };

        if (match.at(i).first == match.at(i).second) {
            i--;
            continue;
        }
        i++;
    }
}

bool CheckTrivially(const CPI& cpi, Match& match, std::map<vertex_t, vertex_t> const& parent,
                    std::set<vertex_t> const& core, std::vector<vertex_t> const& seq) {
    for (std::size_t k = core.size(); k < seq.size(); ++k) {
        std::pair<vertex_t, vertex_t> edge(parent.at(seq.at(k)), seq.at(k));
        std::size_t index = std::find(seq.begin(), seq.end(), parent.at(seq.at(k))) - seq.begin();
        std::pair<std::set<vertex_t>::iterator, std::set<vertex_t>::iterator> its(
                cpi.at(edge).at(*match.at(index).first).begin(),
                cpi.at(edge).at(*match.at(index).first).end());
        match[k] = its;

        while ((match.at(k).first != match.at(k).second) && Visited(match, *match.at(k).first, k)) {
            match.at(k).first++;
        }
        if (match.at(k).first == match.at(k).second) {
            // std::cout << "Trivially satisfied" << std::endl;
            return true;
        }
    }
    return false;
}

bool CheckMatch(const CPI& cpi, Match& match, std::map<vertex_t, vertex_t> const& parent,
                std::set<vertex_t> const& core, std::vector<vertex_t> const& seq,
                graph_t const& graph, graph_t const& query, Gfd const& gfd, int& amount) {
    while (true) {
        std::size_t j = seq.size() - 1;
        while ((j != seq.size()) && (j != core.size() - 1)) {
            if (match.at(j).first == match.at(j).second) {
                std::pair<vertex_t, vertex_t> edge(parent.at(seq.at(j)), seq.at(j));
                std::size_t index =
                        std::find(seq.begin(), seq.end(), parent.at(seq.at(j))) - seq.begin();
                std::pair<std::set<vertex_t>::iterator, std::set<vertex_t>::iterator> its(
                        cpi.at(edge).at(*match.at(index).first).begin(),
                        cpi.at(edge).at(*match.at(index).first).end());
                match[j] = its;
            } else {
                match.at(j).first++;
            }

            while ((match.at(j).first != match.at(j).second) &&
                   Visited(match, *match.at(j).first, j)) {
                match.at(j).first++;
            };
            if (match.at(j).first == match.at(j).second) {
                j--;
                continue;
            }
            j++;
        }
        if (j == core.size() - 1) {
            break;
        }

        amount++;
        // check
        if (!Satisfied(graph, query, seq, match, gfd.GetPremises())) {
            continue;
        }
        if (!Satisfied(graph, query, seq, match, gfd.GetConclusion())) {
            // std::cout << "Checked embeddings: " << amount << std::endl;
            return false;
        }
    }
    return true;
}

bool Check(CPI& cpi, graph_t const& graph, Gfd const& gfd, std::set<vertex_t> const& core,
           std::vector<std::set<vertex_t>> const& forest,
           std::map<vertex_t, vertex_t> const& parent, std::set<edge_t> const& nte) {
    graph_t query = gfd.GetPattern();
    std::vector<std::vector<vertex_t>> paths = GetPaths(core, parent);

    std::vector<vertex_t> nts = {};
    FullNTs(paths, nte, query, nts);
    std::vector<vertex_t> seq = MatchingOrder(cpi, paths, nts);

    CompleteSeq(cpi, forest, parent, query, nte, seq);

    std::set<vertex_t> root_candidates = {};
    std::pair<vertex_t, vertex_t> edge(*seq.begin(), *(++seq.begin()));
    for (auto& vertices : cpi.at(edge)) {
        root_candidates.insert(vertices.first);
    }

    std::vector<std::pair<std::set<vertex_t>::iterator, std::set<vertex_t>::iterator>> match = {};

    if (FullMatch(cpi, match, root_candidates, core, seq, parent, graph, query)) {
        return true;
    }
    int amount = 1;
    // check
    if (Satisfied(graph, query, seq, match, gfd.GetPremises()) &&
        !Satisfied(graph, query, seq, match, gfd.GetConclusion())) {
        // std::cout << "Checked embeddings: " << amount << std::endl;
        return false;
    }

    while (true) {
        int i = static_cast<int>(core.size()) - 1;
        IncrementMatch(i, cpi, match, parent, core, seq, graph, query);
        if (i == -1) {
            break;
        }

        if (forest.empty()) {
            amount++;
            // check
            if (!Satisfied(graph, query, seq, match, gfd.GetPremises())) {
                continue;
            }
            if (!Satisfied(graph, query, seq, match, gfd.GetConclusion())) {
                // std::cout << "Checked embeddings: " << amount << std::endl;
                return false;
            }
            continue;
        }

        if (CheckTrivially(cpi, match, parent, core, seq)) {
            return true;
        }

        if (!CheckMatch(cpi, match, parent, core, seq, graph, query, gfd, amount)) {
            return false;
        }
    }
    // std::cout << "total number of embeddings: " << amount << std::endl;
    return true;
}

bool Validate(graph_t const& graph, Gfd const& gfd) {
    auto start_time = std::chrono::system_clock::now();

    graph_t pat = gfd.GetPattern();
    std::set<std::string> graph_labels = {};
    std::set<std::string> pat_labels = {};
    typename boost::graph_traits<graph_t>::vertex_iterator it, end;
    for (boost::tie(it, end) = vertices(graph); it != end; ++it) {
        graph_labels.insert(graph[*it].attributes.at("label"));
    }
    for (boost::tie(it, end) = vertices(pat); it != end; ++it) {
        pat_labels.insert(pat[*it].attributes.at("label"));
    }
    for (auto const& label : pat_labels) {
        if (graph_labels.find(label) == graph_labels.end()) {
            return true;
        }
    }

    std::set<vertex_t> core = {};
    std::vector<std::set<vertex_t>> forest = {};
    CfDecompose(pat, core, forest);

    int root = GetRoot(graph, pat, core);
    std::vector<std::set<vertex_t>> levels = {};
    std::map<vertex_t, vertex_t> parent;
    std::set<edge_t> snte = {};
    std::set<edge_t> nte = {};
    BfsTree(pat, root, levels, parent, nte, snte);

    std::map<vertex_t, std::set<vertex_t>> candidates;
    CPI cpi;
    TopDownConstruct(cpi, graph, pat, levels, parent, candidates, snte);
    BottomUpRefinement(cpi, graph, pat, levels, parent, candidates);
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);

    // std::cout << "CPI constructed in " << elapsed_milliseconds.count() << ". Matching..."
    //           << std::endl;
    return Check(cpi, graph, gfd, core, forest, parent, nte);
}

}  // namespace

namespace algos {

std::vector<Gfd> EGfdValidation::GenerateSatisfiedGfds(graph_t const& graph,
                                                       std::vector<Gfd> const& gfds) {
    for (auto& gfd : gfds) {
        if (Validate(graph, gfd)) {
            result_.push_back(gfd);
        }
    }
    return result_;
}

}  // namespace algos
