#include "bind_fd.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/fd/fd.h"
#include "algorithms/fd/fd_algorithm.h"
#include "algorithms/fd/mining_algorithms.h"
#include "config/indices/type.h"
#include "py_util/bind_primitive.h"
#include "util/bitset_utils.h"

namespace {
namespace py = pybind11;

template <typename ElementType>
py::tuple VectorToTuple(std::vector<ElementType> vec) {
    std::size_t const size = vec.size();
    py::tuple tuple(size);
    for (std::size_t i = 0; i < size; ++i) {
        tuple[i] = std::move(vec[i]);
    }
    return tuple;
}

py::tuple MakeFdNameTuple(FD const& fd) {
    auto [lhs, rhs] = fd.ToNameTuple();
    return py::make_tuple(VectorToTuple(std::move(lhs)), std::move(rhs));
}
}  // namespace

namespace python_bindings {
void BindFd(py::module_& main_module) {
    using namespace algos;

    auto fd_module = main_module.def_submodule("fd");
    py::class_<FD>(fd_module, "FD")
            .def("__str__", &FD::ToLongString)
            .def("to_long_string", &FD::ToLongString)
            .def("to_short_string", &FD::ToShortString)
            .def("to_index_tuple",
                 [](FD const& fd) {
                     return py::make_tuple(VectorToTuple(fd.GetLhsIndices()),
                                           std::move(fd.GetRhsIndex()));
                 })
            .def("to_name_tuple", MakeFdNameTuple)
            // TODO: implement proper equality check for FD
            .def("__eq__", [](FD const& fd1,
                              FD const& fd2) { return fd1.ToNameTuple() == fd2.ToNameTuple(); })
            .def("__hash__", [](FD const& fd) { return py::hash(MakeFdNameTuple(fd)); })
            .def_property_readonly("lhs_indices", &FD::GetLhsIndices)
            .def_property_readonly("rhs_index", &FD::GetRhsIndex);

    static constexpr auto kPyroName = "Pyro";
    static constexpr auto kTaneName = "Tane";
    static constexpr auto kPFDTaneName = "PFDTane";
    auto fd_algos_module =
            BindPrimitive<hyfd::HyFD, Aid, Depminer, DFD, FastFDs, FDep, FdMine, FUN, Pyro, Tane,
                          PFDTane>(fd_module, &FDAlgorithm::SortedFdList, "FdAlgorithm", "get_fds",
                                   {"HyFD", "Aid", "Depminer", "DFD", "FastFDs", "FDep", "FdMine",
                                    "FUN", kPyroName, kTaneName, kPFDTaneName});

    auto define_submodule = [&fd_algos_module, &main_module](char const* name,
                                                             std::vector<char const*> algorithms) {
        auto algos_module = main_module.def_submodule(name).def_submodule("algorithms");
        for (auto algo_name : algorithms) {
            algos_module.attr(algo_name) = fd_algos_module.attr(algo_name);
        }
        algos_module.attr("Default") = algos_module.attr(algorithms.front());
    };

    define_submodule("afd", {kPyroName, kTaneName});
    define_submodule("pfd", {kPFDTaneName});
}
}  // namespace python_bindings
