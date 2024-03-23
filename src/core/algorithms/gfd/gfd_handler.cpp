#include "gfd_handler.h"

#include <iostream>
#include <set>
#include <sstream>
#include <thread>

#include <boost/graph/eccentricity.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "balancer.h"
#include "config/equal_nulls/option.h"
#include "config/names_and_descriptions.h"
#include "config/option_using.h"
#include "config/tabular_data/input_table/option.h"
#include "config/thread_number/option.h"

namespace algos {

std::vector<std::string> GfdHandler::StringGfdList() {
    std::vector<std::string> result;
    result.reserve(result_.size());
    for (auto& gfd : result_) {
        std::stringstream gfd_stream;
        parser::graph_parser::WriteGfd(gfd_stream, gfd);
        result.push_back(gfd_stream.str());
    }
    return result;
}

GfdHandler::GfdHandler() : Algorithm({}) {
    RegisterOptions();
    MakeOptionsAvailable({config::names::kGfdData, config::names::kGraphData});
};

void GfdHandler::RegisterOptions() {
    using namespace config::names;
    using namespace config::descriptions;
    DESBORDANTE_OPTION_USING;

    RegisterOption(config::Option{&gfd_paths_, kGfdData, kDGfdData});
    RegisterOption(config::Option{&graph_path_, kGraphData, kDGraphData});
}

void GfdHandler::LoadDataInternal() {
    std::ifstream f(graph_path_);
    graph_ = parser::graph_parser::ReadGraph(f);
    f.close();
    for (auto const& path : gfd_paths_) {
        auto gfd_path = path;
        f.open(gfd_path);
        Gfd gfd = parser::graph_parser::ReadGfd(f);
        f.close();
        gfds_.push_back(gfd);
    }
}

void GfdHandler::ResetState() {}

unsigned long long GfdHandler::ExecuteInternal() {
    auto start_time = std::chrono::system_clock::now();

    result_ = GenerateSatisfiedGfds(graph_, gfds_);

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    // std::cout << "Satisfied GFDs: " << result_.size() << "/" << gfds_.size() << std::endl;
    return elapsed_milliseconds.count();
}

}  // namespace algos
