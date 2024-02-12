#include "bind_cfd.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/cfd/fd_first_algorithm.h"
#include "algorithms/cfd/model/raw_cfd.h"
#include "py_util/bind_primitive.h"

namespace py = pybind11;

namespace python_bindings {
void BindCfd(py::module_& main_module) {
    using namespace algos::cfd;

    auto cfd_module = main_module.def_submodule("cfd");

    py::class_<RawCFD::RawItem>(cfd_module, "Item")
            .def_property_readonly("attribute", &RawCFD::RawItem::GetAttribute)
            .def_property_readonly("value", &RawCFD::RawItem::GetValue);

    py::class_<RawCFD>(cfd_module, "CFD")
            .def("__str__", &RawCFD::ToString)
            .def_property_readonly("lhs_items", &RawCFD::GetLhs)
            .def_property_readonly("rhs_item", &RawCFD::GetRhs);

    BindPrimitive<FDFirstAlgorithm>(cfd_module, &CFDDiscovery::GetCfds, "CfdAlgorithm", "get_cfds",
                                    {"FDFirst"});
}
}  // namespace python_bindings
