#pragma once
#include <vector>

#include "algorithms/algorithm.h"
#include "config/names_and_descriptions.h"
#include "gfd.h"
#include "parser/graph_parser/graph_parser.h"

namespace algos {

class GfdHandler : public Algorithm {
protected:
    std::filesystem::path graph_path_;
    std::vector<std::filesystem::path> gfd_paths_;

    graph_t graph_;
    std::vector<Gfd> gfds_;
    std::vector<Gfd> result_;

    unsigned long long ExecuteInternal();

    void ResetState() final;
    void LoadDataInternal() final;

    void RegisterOptions();

public:
    virtual std::vector<Gfd> GenerateSatisfiedGfds(graph_t const& graph,
                                                   std::vector<Gfd> const& gfds) = 0;

    GfdHandler();

    GfdHandler(graph_t graph_, std::vector<Gfd> gfds_)
        : Algorithm({}), graph_(graph_), gfds_(gfds_) {
        ExecutePrepare();
    }

    std::vector<Gfd> GfdList() {
        return result_;
    }

    std::vector<std::string> StringGfdList();
};

}  // namespace algos
