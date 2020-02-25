#include "EmpsvdCore.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <cstddef>

namespace py = pybind11;


PYBIND11_MODULE(pyempsvd, mod) {
	mod.doc() = "EM-PSVD module implimented by C++11";

	py::class_<Empsvd::EmpsvdCore>(mod, "EmpsvdCore")
	  .def(py::init<size_t, const Eigen::ArrayXd&, const Eigen::ArrayXd&, const Eigen::ArrayXXd&, size_t, double, bool, bool>(),
			py::arg("k"), py::arg("x"), py::arg("y"), py::arg("theta0"),
			py::arg("max_iter") = 1000, py::arg("tol") = 1e-2, py::arg("fix_alpha") = false, py::arg("fix_ab") = false)
		.def(py::init<size_t, const Eigen::ArrayXd&, const Eigen::ArrayXd&, size_t, double, bool, bool>(),
			py::arg("k"), py::arg("x"), py::arg("y"),
			py::arg("max_iter") = 1000, py::arg("tol") = 1e-2, py::arg("fix_alpha") = false, py::arg("fix_ab") = false)
	  .def(py::init<size_t, const Eigen::ArrayXd&, const Eigen::ArrayXd&, const Eigen::ArrayXd&, const Eigen::ArrayXXd&, size_t, double, bool, bool>(),
			py::arg("k"), py::arg("x"), py::arg("y"), py::arg("z"), py::arg("theta0"),
			py::arg("max_iter") = 1000, py::arg("tol") = 1e-2, py::arg("fix_alpha") = false, py::arg("fix_ab") = false)
		.def(py::init<size_t, const Eigen::ArrayXd&, const Eigen::ArrayXd&, const Eigen::ArrayXd&, size_t, double, bool, bool>(),
			py::arg("k"), py::arg("x"), py::arg("y"), py::arg("z"),
			py::arg("max_iter") = 1000, py::arg("tol") = 1e-2, py::arg("fix_alpha") = false, py::arg("fix_ab") = false)
		.def("fit", &Empsvd::EmpsvdCore::fit)
		.def("e_step", &Empsvd::EmpsvdCore::e_step)
		.def("m_step", &Empsvd::EmpsvdCore::m_step)
		.def("get_aic", &Empsvd::EmpsvdCore::get_aic)
		.def("get_bic", &Empsvd::EmpsvdCore::get_bic)
		.def("get_loglikelihood", (double (Empsvd::EmpsvdCore::*)(void)) & Empsvd::EmpsvdCore::get_loglikelihood)
		.def("get_loglikelihood", (double (Empsvd::EmpsvdCore::*)(const Eigen::ArrayXXd&)) &Empsvd::EmpsvdCore::get_loglikelihood, py::arg("theta"))

		.def_readonly("x", &Empsvd::EmpsvdCore::x)
		.def_readonly("y", &Empsvd::EmpsvdCore::y)
		.def_readonly("z", &Empsvd::EmpsvdCore::z)
		.def_readwrite("theta", &Empsvd::EmpsvdCore::theta)
		.def_readwrite("gamma", &Empsvd::EmpsvdCore::gamma)
		.def_readonly("k", &Empsvd::EmpsvdCore::k)
		.def_readonly("n", &Empsvd::EmpsvdCore::n)
		.def_property_readonly_static("m", [](py::object) {return Empsvd::EmpsvdCore::m; })
		.def_readwrite("niter", &Empsvd::EmpsvdCore::niter)

		.def_static("make_theta0", (Eigen::ArrayXXd (*)(const size_t, const Eigen::ArrayXd&, const Eigen::ArrayXd&, const Eigen::ArrayXd&)) &Empsvd::EmpsvdCore::make_theta0, py::arg("k"), py::arg("x"), py::arg("y"), py::arg("z"))
		.def_static("make_theta0", (Eigen::ArrayXXd (*)(const size_t, const Eigen::ArrayXd&, const Eigen::ArrayXd&)) &Empsvd::EmpsvdCore::make_theta0, py::arg("k"), py::arg("x"), py::arg("y"))
		.def_static("calc_log_pxy", &Empsvd::EmpsvdCore::calc_log_pxy, py::arg("x"), py::arg("y"), py::arg("theta"))
		.def_static("calc_pxy", &Empsvd::EmpsvdCore::calc_pxy, py::arg("x"), py::arg("y"), py::arg("theta"))
		;

}
