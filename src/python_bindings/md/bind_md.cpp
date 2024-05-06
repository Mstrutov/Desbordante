#include "md/bind_md.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/md/hymd/preprocessing/similarity_measure/levenshtein_similarity_measure.h"
#include "algorithms/md/hymd/similarity_measure_creator.h"
#include "algorithms/md/md.h"
#include "algorithms/md/mining_algorithms.h"
#include "py_util/bind_primitive.h"

namespace {
namespace py = pybind11;
using model::MD;
}  // namespace

namespace python_bindings {
void BindMd(py::module_& main_module) {
    using namespace algos;
    using namespace algos::hymd;
    using namespace algos::hymd::preprocessing::similarity_measure;
    using namespace py::literals;

    auto md_module = main_module.def_submodule("md");
    py::class_<MD>(md_module, "MD")
            .def_property_readonly("lhs_bounds", &MD::GetLhsDecisionBounds)
            .def_property_readonly("rhs", &MD::GetRhs)
            .def("to_long_string", &MD::ToStringFull)
            .def("to_short_string", &MD::ToStringShort)
            .def("__str__", &MD::ToStringShort);
    auto measures_module = md_module.def_submodule("similarity_measures");
    py::class_<SimilarityMeasureCreator, std::shared_ptr<SimilarityMeasureCreator>>(
            measures_module, "SimilarityMeasure");
    py::class_<LevenshteinSimilarityMeasure::Creator, SimilarityMeasureCreator,
               std::shared_ptr<LevenshteinSimilarityMeasure::Creator>>(measures_module,
                                                                       "LevenshteinSimilarity")
            .def(py::init<model::md::DecisionBoundary, bool, std::size_t>(),
                 "minimum_similarity"_a = 0.7, "is_null_equal_null"_a = true,
                 "bound_number_limit"_a = 0);
    BindPrimitive<HyMD>(md_module, &MdAlgorithm::MdList, "MdAlgorithm", "get_mds", {"HyMD"});
}
}  // namespace python_bindings
