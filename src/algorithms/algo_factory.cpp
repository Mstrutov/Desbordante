#include "algo_factory.h"

namespace algos {

std::unique_ptr<FDAlgorithm> CreateFDAlgorithmInstance(Algo const algo,
                                                       std::string const& dataset,
                                                       char const sep,
                                                       bool const has_header,
                                                       double const error,
                                                       unsigned int const max_lhs,
                                                       unsigned int const parallelism,
                                                       int seed) {
    auto path = std::filesystem::current_path() / "inputData" / dataset;

    /* A little bit ugly, but I don't know how to do it better */
    switch (algo) {
    case Algo::depminer:
        return std::make_unique<Depminer>(path, sep, has_header);
    case Algo::dfd:
        return std::make_unique<DFD>(path, sep, has_header, parallelism);
    case Algo::fastfds:
        return std::make_unique<FastFDs>(path, sep, has_header, max_lhs, parallelism);
    case Algo::fdep:
        return std::make_unique<FDep>(path, sep, has_header);
    case Algo::fdmine:
        return std::make_unique<Fd_mine>(path, sep, has_header);
    case Algo::pyro:
        return std::make_unique<Pyro>(path, sep, has_header, seed, error, max_lhs, parallelism);
    case Algo::tane:
        return std::make_unique<Tane>(path, sep, has_header, error, max_lhs);
    default:
        assert(0);
        break;
    }

    return nullptr;
}

} // namespace algos

