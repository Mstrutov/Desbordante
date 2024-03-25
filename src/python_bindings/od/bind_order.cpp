#include "bind_order.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "algorithms/od/order/list_od.h"
#include "algorithms/od/order/order.h"
#include "py_util/bind_primitive.h"

namespace {
namespace py = pybind11;
}  // namespace

namespace python_bindings {

void BindOrder(py::module_& main_module) {
    using namespace algos::order;

    auto od_module = main_module.def_submodule("od");

    py::class_<ListOD>(od_module, "ListOD")
            .def_readonly("lhs", &ListOD::lhs)
            .def_readonly("rhs", &ListOD::rhs);

    static constexpr auto kOrderName = "Order";

    auto od_algos_module =
            BindPrimitiveNoBase<Order>(od_module, "Order").def("get_list_ods", [](Order& algo) {
                OrderDependencies const& map_res = algo.GetValidODs();
                std::vector<ListOD> res;
                for (auto const& [lhs, rhs_list] : map_res) {
                    for (AttributeList const& rhs : rhs_list) {
                        res.emplace_back(lhs, rhs);
                    }
                }
                return res;
            });

    main_module.attr("od_module") = od_module;
}

}  // namespace python_bindings
