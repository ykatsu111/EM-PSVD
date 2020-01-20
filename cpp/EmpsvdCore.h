#pragma once

#include <Eigen/Dense>
#include <cstddef>

namespace Empsvd {

	class EmpsvdCore
	{
	public:
		const Eigen::ArrayXd x, y;
		Eigen::ArrayXXd theta;
		Eigen::ArrayXXd gamma;
		const size_t k, n;
		static const size_t m = 6;
		const size_t max_iter;
		size_t niter;
		const double tol;
		const bool fix_alpha, fix_ab;

		// もしかして、参照渡しは危険？ユーザーが元のデータを消してるけどこのクラスオブジェクトが残っている状況がありそう？
		EmpsvdCore(
			Eigen::ArrayXd x, Eigen::ArrayXd y, size_t k, Eigen::ArrayXXd theta0,
			size_t max_iter = 1000, double tol = 1e-2, bool fix_alpha = false, bool fix_ab = false
		);
		EmpsvdCore(
			Eigen::ArrayXd x, Eigen::ArrayXd y, size_t k,
			size_t max_iter = 1000, double tol = 1e-2, bool fix_alpha = false, bool fix_ab = false
		);
		~EmpsvdCore();

		void fit();
		void e_step();
		void m_step();

		double get_aic();
		double get_bic();
		double get_loglikelihood();
		double get_loglikelihood(const Eigen::ArrayXXd& theta);

		static Eigen::ArrayXXd make_theta0(
			const Eigen::ArrayXd& x, const Eigen::ArrayXd& y,
			size_t const k, size_t const m
		);
		static Eigen::ArrayXXd calc_log_pxy(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXXd& theta);
		static Eigen::ArrayXXd calc_pxy(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXXd& theta);
		static double digammad(double a);

	private:
		void check_init();
		Eigen::ArrayXXd get_log_pxy();
		Eigen::ArrayXXd get_log_pxy(const Eigen::ArrayXXd& theta);
		Eigen::ArrayXXd get_pxy();
		Eigen::ArrayXXd get_pxy(const Eigen::ArrayXXd& theta);
		Eigen::ArrayXd get_logsum_pxy();
		Eigen::ArrayXd get_logsum_pxy(Eigen::ArrayXXd log_pxy);
		Eigen::ArrayXd get_sum_pxy();
		Eigen::ArrayXd get_sum_pxy(Eigen::ArrayXXd log_pxy);
		Eigen::ArrayXXd get_gamma();
		Eigen::ArrayXXd get_gamma(const Eigen::ArrayXXd& theta);
		double get_new_omegak(Eigen::Index ik);
		double get_new_ak(Eigen::Index ik, double new_bk);
		double get_new_bk(Eigen::Index ik);
		double get_new_sigma2k(Eigen::Index ik, double new_ak, double new_bk);
		double get_new_alphak(Eigen::Index ik);
		double get_new_lambdak(Eigen::Index ik, double new_alphak);
		double bkdot(Eigen::Index ik, double bk);
		double get_new_bk_by_newton(Eigen::Index ik, double bk0, double bk1, size_t max_iter = 100, double tol = 1e-8);
		double get_new_alphak_by_invdigamma(Eigen::Index ik, double alphak0, size_t max_iter = 500, double tol = 1e-2);
	};

}
