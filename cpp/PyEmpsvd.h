#pragma once

#include "EmpsvdCore.h"
#include <pybind11/eigen.h>
#include <cstddef>

namespace Empsvd {

	class PyEmpsvd : public Empsvd::EmpsvdCore
	{
	public:
		const Eigen::ArrayXd& getX();
		const Eigen::ArrayXd& getY();
		Eigen::ArrayXXd& getTheta();
		Eigen::ArrayXXd& getGamma();
		size_t getK();
		size_t getN();
		size_t getM();
		const size_t& getNiter();
	};

}