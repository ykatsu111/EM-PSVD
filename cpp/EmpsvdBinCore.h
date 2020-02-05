#pragma once

#include <Eigen/Dense>
#include <cstddef>
#include "EmpsvdCore.h"

namespace Empsvd {

	class EmpsvdBinCore : public EmpsvdCore
	{
	public:
		const Eigen::ArrayXd z;

		EmpsvdBinCore(
			Eigen::ArrayXd x, Eigen::ArrayXd y, Eigen::ArrayXd z, size_t k, Eigen::ArrayXXd theta0,
			size_t max_iter = 1000, double tol = 1e-2, bool fix_alpha = false, bool fix_ab = false
		);
		EmpsvdBinCore(
			Eigen::ArrayXd x, Eigen::ArrayXd y, Eigen::ArrayXd z, size_t k,
			size_t max_iter = 1000, double tol = 1e-2, bool fix_alpha = false, bool fix_ab = false
		);
		
		using EmpsvdCore::get_loglikelihood;
		double get_loglikelihood(const Eigen::ArrayXXd& theta) override;

		static Eigen::ArrayXXd make_theta0(
			const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z,
			size_t const k, size_t const m
		);		

	protected:
		using EmpsvdCore::get_gamma;
		Eigen::ArrayXXd get_gamma(const Eigen::ArrayXXd& theta) override;

	};

}
