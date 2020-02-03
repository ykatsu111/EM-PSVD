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
			Eigen::ArrayXXd x, Eigen::ArrayXXd y, Eigen::ArrayXXd z, size_t k, Eigen::ArrayXXd theta0,
			size_t max_iter = 1000, double tol = 1e-2, bool fix_alpha = false, bool fix_ab = false
		);
		EmpsvdBinCore(
			Eigen::ArrayXXd x, Eigen::ArrayXXd y, Eigen::ArrayXXd z, size_t k,
			size_t max_iter = 1000, double tol = 1e-2, bool fix_alpha = false, bool fix_ab = false
		);

		static Eigen::ArrayXXd make_theta0(
			const Eigen::ArrayXXd& x, const Eigen::ArrayXXd& y, const Eigen::ArrayXXd& z,
			size_t const k, size_t const m
		);
		static Eigen::ArrayXXd calc_log_pxy(
			const Eigen::ArrayXXd& x, const Eigen::ArrayXXd& y, const Eigen::ArrayXXd& z,
			const Eigen::ArrayXXd& theta
		);
		static Eigen::ArrayXXd calc_pxy(
			const Eigen::ArrayXXd& x, const Eigen::ArrayXXd& y, const Eigen::ArrayXXd& z,
			const Eigen::ArrayXXd& theta
		);

	private:
		Eigen::ArrayXXd get_log_pxy() override;
		Eigen::ArrayXXd Empsvd::EmpsvdCore::get_log_pxy(const Eigen::ArrayXXd& theta) override;
		double get_new_omegak(Eigen::Index ik) override;
		double get_new_ak(Eigen::Index ik, double new_bk) override;
		double get_new_sigma2k(Eigen::Index ik, double new_ak, double new_bk) override;
		double get_new_lambdak(Eigen::Index ik, double new_alphak) override;
		double bkdot(Eigen::Index ik, double bk) override;
		double get_new_alphak_by_invdigamma(Eigen::Index ik, double alphak0, size_t max_iter = 500, double tol = 1e-2) override;
	};

}