#include "EmpsvdBinCore.h"
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <cstddef>

Empsvd::EmpsvdBinCore::EmpsvdBinCore(
	Eigen::ArrayXd x, Eigen::ArrayXd y, Eigen::ArrayXd z, size_t k, Eigen::ArrayXXd theta0, 
	size_t max_iter, double tol, bool fix_alpha, bool fix_ab
) : Empsvd::EmpsvdCore(x, y, k, theta0,	max_iter, tol, fix_alpha, fix_ab), z(z)
{
}

Empsvd::EmpsvdBinCore::EmpsvdBinCore(
	Eigen::ArrayXd x, Eigen::ArrayXd y, Eigen::ArrayXd z, size_t k, 
	size_t max_iter, double tol, bool fix_alpha, bool fix_ab
) : Empsvd::EmpsvdCore(x, y, k, this->make_theta0(x, y, z, k, this->m),	max_iter, tol, fix_alpha, fix_ab), z(z)
{
}

Eigen::ArrayXXd Empsvd::EmpsvdBinCore::make_theta0(
	const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z, 
	size_t const k, size_t const m
)
{
	Eigen::ArrayXXd theta0(k, m);

	int const n = z.sum();
	Eigen::ArrayXd p = z / n;
	double const y_mean = (p * y).sum();
	double const y_var = ((y - y_mean).square() * p).sum();
	double const x_mean = (p * x).sum();
	double const x_var = ((x - x_mean).square() * p).sum();

	theta0.col(0) = 1. / k;                        // omega: mixing fraction
	theta0.col(3) = y_var;                         // sigma2: variance of velocity
	theta0.col(4) = std::pow(x_mean, 2) / x_var;   // alpha: mu + 1, shape parameter of size distribution
	theta0.col(5) = x_mean / x_var;                // lambda: slope parameter of size distribution

	Eigen::ArrayXd exist_data = (z.array() > 0.).select(
		Eigen::ArrayXd::Constant(z.size(), 1),
		Eigen::ArrayXd::Constant(z.size(), 0)
	);
	double a_upp = (y * exist_data).maxCoeff() / x_mean;
	double a_low = y_mean / 4;

	if (k > 1) {
		Eigen::ArrayXd lins = Eigen::ArrayXd::LinSpaced(k, 0., 1.);
		theta0.col(1) = (lins * std::log(a_upp / a_low) + std::log(a_low)).exp();  // a: an intercept parameter of v-d relationship
		theta0.col(2) = lins;                                                      // b: a slope parameter of v-d relationship
	}
	else {
		theta0.col(1) = std::exp(std::log(a_upp * a_low) / 2.);
		theta0.col(2) = 0.5;
	}

	return theta0;
}

Eigen::ArrayXXd Empsvd::EmpsvdBinCore::calc_log_pxy(
	const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z, const Eigen::ArrayXXd& theta
)
{
	Eigen::ArrayXXd log_pxy = calc_log_pxy(x, y, theta);
	for (Eigen::Index ik = 0; ik < theta.rows(); ik++) {
		log_pxy *= z;
	}
	return log_pxy;
}

Eigen::ArrayXXd Empsvd::EmpsvdBinCore::calc_pxy(
	const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z, const Eigen::ArrayXXd& theta
)
{
	return calc_log_pxy(x, y, z, theta).exp();
}

Eigen::ArrayXXd Empsvd::EmpsvdBinCore::get_log_pxy()
{
	return this->get_log_pxy(this->theta);
}

Eigen::ArrayXXd Empsvd::EmpsvdCore::get_log_pxy(const Eigen::ArrayXXd& theta)
{
	return this->calc_log_pxy(this->x, this->y, this->z, theta);
}