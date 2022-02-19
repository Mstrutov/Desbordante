#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <easylogging++.h>

#include "AlgoFactory.h"

namespace po = boost::program_options;

INITIALIZE_EASYLOGGINGPP

template<typename BetterEnumType>
static std::string EnumToAvailableValues() {
    std::stringstream avail_values;

    avail_values << '[';

    for (auto const& name : BetterEnumType::_names()) {
        avail_values << name << '|';
    }

    avail_values.seekp(-1, avail_values.cur);
    avail_values << ']';

    return avail_values.str();
}

static bool CheckOptions(std::string const& task, std::string const& alg, double error) {
    if (!algos::AlgoMiningType::_is_valid(task.c_str())) {
        std::cout << "ERROR: no matching task."
                     " Available tasks (primitives to mine) are:\n" +
                     EnumToAvailableValues<algos::AlgoMiningType>() + '\n';
        return false;
    }

    if (!algos::Algo::_is_valid(alg.c_str())) {
        std::cout << "ERROR: no matching algorithm."
                     " Available algorithms are:\n" +
                     EnumToAvailableValues<algos::Algo>() + '\n';
        return false;
    }

    if (error > 1 || error < 0) {
        std::cout << "ERROR: error should be between 0 and 1.\n";
    }

    return true;
}

int main(int argc, char const* argv[]) {
    std::string algo;
    std::string dataset;
    std::string task;
    char separator = ',';
    bool has_header = true;
    int seed = 0;
    double error = 0.0;
    unsigned int max_lhs = -1;
    ushort threads = 0;
    bool is_null_equal_null = true;

    std::string const algo_desc = "algorithm to use. Available algorithms:\n" +
                                  EnumToAvailableValues<algos::Algo>() +
                                  " for FD mining.";
    std::string const task_desc = "type of dependency to mine. Available tasks:\n" +
                                  EnumToAvailableValues<algos::AlgoMiningType>();

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "print this message and exit")
        ("task", po::value<std::string>(&task),
         task_desc.c_str())
        ("algo", po::value<std::string>(&algo),
         algo_desc.c_str())
        ("data", po::value<std::string>(&dataset),
         "path to CSV file, relative to ./inputData")
        ("separator,s", po::value<char>(&separator)->default_value(separator),
         "CSV separator")
        ("has_header", po::value<bool>(&has_header)->default_value(has_header),
         "CSV header presence flag [true|false]")
        ("seed", po::value<int>(&seed)->default_value(seed), "RNG seed")
        ("error", po::value<double>(&error)->default_value(error),
         "error for AFD algorithms")
        ("max_lhs", po::value<unsigned int>(&max_lhs)->default_value(max_lhs),
         "max considered LHS size")
        ("threads", po::value<ushort>(&threads)->default_value(threads),
         "number of threads to use. If 0 is specified then as many threads are used as "
         "the hardware can handle concurrently.")
         ("is_null_equal_null", po::value<bool>(&is_null_equal_null)->default_value(true),
          "Is NULL value equals another NULL value")
        ;

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (po::error &e) {
        std::cout << e.what() << std::endl;
        return 0;
    }

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 0;
    }

    el::Loggers::configureFromGlobal("logging.conf");

    std::transform(algo.begin(), algo.end(), algo.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (!CheckOptions(task, algo, error)) {
        std::cout << desc << std::endl;
        return 1;
    }

    /* Remove options that are not related to the algorithm configuration */
    vm.erase("task");
    vm.erase("algo");

    std::cout << "Input: algorithm \"" << algo
              << "\" with seed " << std::to_string(seed)
              << ", error \"" << std::to_string(error)
              << ", max_lhs \"" << std::to_string(max_lhs)
              << "\" and dataset \"" << dataset
              << "\" with separator \'" << separator
              << "\'. Header is " << (has_header ? "" : "not ") << "present. " << std::endl;

    std::unique_ptr<algos::Primitive> algorithm_instance =
        algos::CreateAlgorithmInstance(task, algo, vm);

    try {
        unsigned long long elapsed_time = algorithm_instance->Execute();
        std::cout << "> ELAPSED TIME: " << elapsed_time << std::endl;
    } catch (std::runtime_error& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
