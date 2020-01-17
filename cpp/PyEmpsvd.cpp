#include "PyEmpsvd.h"
#include "EmpsvdCore.h"
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <cstddef>

namespace py = pybind11;


PYBIND11_MODULE(pyempsvd, mod) {
	mod.doc() = "EM-PSVD module implimented by C++11";

	py::class_<Empsvd::EmpsvdCore>(mod, "Empsvd")
		.def(py::init<Eigen::ArrayXd, Eigen::ArrayXd, size_t, Eigen::ArrayXd, size_t, double, bool, bool>(),
			py::arg("x"), py::arg("y"), py::arg("k"), py::arg("theta0"),
			py::arg("max_iter") = 1000, py::arg("tol") = 1e-2, py::arg("fix_alpha") = false, py::arg("fix_ab") = false)
		.def(py::init<Eigen::ArrayXd, Eigen::ArrayXd, size_t, size_t, double, bool, bool>(),
			py::arg("x"), py::arg("y"), py::arg("k"),
			py::arg("max_iter") = 1000, py::arg("tol") = 1e-2, py::arg("fix_alpha") = false, py::arg("fix_ab") = false)
		.def("fit", &Empsvd::EmpsvdCore::fit)
		.def("e_step", &Empsvd::EmpsvdCore::e_step)
		.def("m_step", &Empsvd::EmpsvdCore::m_step)
		.def("get_aic", &Empsvd::EmpsvdCore::get_aic)
		.def("get_bic", &Empsvd::EmpsvdCore::get_bic)
		.def("get_loglikelihood", (double (Empsvd::EmpsvdCore::*)(void)) & Empsvd::EmpsvdCore::get_loglikelihood)
		.def("get_loglikelihood", (double (Empsvd::EmpsvdCore::*)(const Eigen::ArrayXXd&)) & Empsvd::EmpsvdCore::get_loglikelihood, py::arg("theta"))

		.def_readonly("x", &Empsvd::EmpsvdCore::x)
		.def_readonly("y", &Empsvd::EmpsvdCore::y)
		.def_readwrite("theta", &Empsvd::EmpsvdCore::theta)
		.def_readwrite("gamma", &Empsvd::EmpsvdCore::gamma)
		.def_readonly("k", &Empsvd::EmpsvdCore::k)
		.def_readonly("n", &Empsvd::EmpsvdCore::n)
		.def_property_readonly_static("m", [](py::object) {return Empsvd::EmpsvdCore::m; })
		.def_readwrite("niter", &Empsvd::EmpsvdCore::niter)
		;
}

const Eigen::ArrayXd& Empsvd::PyEmpsvd::getX()
{
	return this->x;
}

const Eigen::ArrayXd& Empsvd::PyEmpsvd::getY()
{
	return this->y;
}

Eigen::ArrayXXd& Empsvd::PyEmpsvd::getTheta()
{
	return this->theta;
}

Eigen::ArrayXXd& Empsvd::PyEmpsvd::getGamma()
{
	return this->gamma;
}

size_t Empsvd::PyEmpsvd::getK()
{
	return this->k;
}

size_t Empsvd::PyEmpsvd::getN()
{
	return this->n;
}

size_t Empsvd::PyEmpsvd::getM()
{
	return this->m;
}

const size_t& Empsvd::PyEmpsvd::getNiter()
{
	return this->niter;
}
