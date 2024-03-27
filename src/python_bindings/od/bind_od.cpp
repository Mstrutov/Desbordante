#include "bind_od.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/od/fastod/fastod.h"
#include "algorithms/od/fastod/model/canonical_od.h"
#include "algorithms/od/mining_algorithms.h"
#include "config/indices/type.h"
#include "py_util/bind_primitive.h"
#include "util/bitset_utils.h"

namespace {
namespace py = pybind11;
}  // namespace

namespace python_bindings {

void BindOd(py::module_& main_module) {
    using namespace algos::fastod;

    auto od_module = main_module.def_submodule("od");

    py::class_<AscCanonicalOD>(od_module, "AscCanonicalOD")
            .def("__str__", &AscCanonicalOD::ToString)
            .def("__eq__",
                 [](AscCanonicalOD const& od1, AscCanonicalOD const& od2) {
                     return od1.ToString() == od2.ToString();
                 })
            .def("__hash__",
                 [](AscCanonicalOD const& od) { return py::hash(py::str(od.ToString())); });

    py::class_<DescCanonicalOD>(od_module, "DescCanonicalOD")
            .def("__str__", &DescCanonicalOD::ToString)
            .def("__eq__",
                 [](DescCanonicalOD const& od1, DescCanonicalOD const& od2) {
                     return od1.ToString() == od2.ToString();
                 })
            .def("__hash__",
                 [](DescCanonicalOD const& od) { return py::hash(py::str(od.ToString())); });

    py::class_<SimpleCanonicalOD>(od_module, "SimpleCanonicalOD")
            .def("__str__", &SimpleCanonicalOD::ToString)
            .def("__eq__",
                 [](SimpleCanonicalOD const& od1, SimpleCanonicalOD const& od2) {
                     return od1.ToString() == od2.ToString();
                 })
            .def("__hash__",
                 [](SimpleCanonicalOD const& od) { return py::hash(py::str(od.ToString())); });

    static constexpr auto kFastodName = "Fastod";

    auto od_algos_module = BindPrimitiveNoBase<algos::Fastod>(od_module, "Fastod")
                                   .def("get_asc_ods", &algos::Fastod::GetAscendingDependencies)
                                   .def("get_desc_ods", &algos::Fastod::GetDescendingDependencies)
                                   .def("get_simple_ods", &algos::Fastod::GetSimpleDependencies);

    main_module.attr("od_module") = od_module;
}

}  // namespace python_bindings
