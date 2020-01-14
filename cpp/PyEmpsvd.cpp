#include "PyEmpsvd.h"
#include "EmpsvdCore.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <cstddef>

namespace py = pybind11;

PYBIND11_MODULE(empsvd, mod) {
	mod.doc() = "EM-PSVD module implimented by C++11";

	py::class_<EmpsvdCore>(mod, "Empsvd")
		.def(py::init<const Eigen::ArrayXd&, const Eigen::ArrayXd&, size_t, Eigen::ArrayXd, size_t, double, bool, bool>(),
			py::arg("x"), py::arg("y"), py::arg("k"), py::arg("theta0"), 
			py::arg("max_iter") = 1000, py::arg("tol") = 1e-2, py::arg("fix_alpha") = false, py::arg("fix_ab") = false)
		.def(py::init<const Eigen::ArrayXd&, const Eigen::ArrayXd&, size_t, size_t, double, bool, bool>(),
			py::arg("x"), py::arg("y"), py::arg("k"),
			py::arg("max_iter") = 1000, py::arg("tol") = 1e-2, py::arg("fix_alpha") = false, py::arg("fix_ab") = false)
		.def("fit", &EmpsvdCore::fit)
		.def("e_step", &EmpsvdCore::e_step)
		.def("m_step", &EmpsvdCore::m_step)
		.def("get_aic", &EmpsvdCore::get_aic, py::return_value_policy::reference_internal)
		.def("get_bic", &EmpsvdCore::get_bic, py::return_value_policy::reference_internal)
		.def("get_loglikelihood", &EmpsvdCore::get_loglikelihood, py::return_value_policy::reference_internal)
		.def("get_loglikelihood", &EmpsvdCore::get_loglikelihood, py::arg("theta"), py::return_value_policy::reference_internal)

}