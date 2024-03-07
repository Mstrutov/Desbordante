#include "bind_ucc_verification.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/ucc/verification_algorithms.h"
#include "py_util/bind_primitive.h"

namespace {
namespace py = pybind11;
}  // namespace

namespace python_bindings {
void BindUccVerification(pybind11::module_& main_module) {
    using namespace algos;

    auto ucc_verification_module = main_module.def_submodule("ucc_verification");

    BindPrimitiveNoBase<UCCVerifier>(ucc_verification_module, "UccVerifier")
            .def("ucc_holds", &UCCVerifier::UCCHolds)
            .def("get_num_clusters_violating_ucc", &UCCVerifier::GetNumClustersViolatingUCC)
            .def("get_num_rows_violating_ucc", &UCCVerifier::GetNumRowsViolatingUCC)
            .def("get_clusters_violating_ucc", &UCCVerifier::GetClustersViolatingUCC)
            .def("get_error", &UCCVerifier::GetError);
    main_module.attr("aucc_verification") = ucc_verification_module;
}
}  // namespace python_bindings
