#include "EmpsvdCore.h"
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <cstddef>


const size_t Empsvd::EmpsvdCore::m;


Empsvd::EmpsvdCore::EmpsvdCore(
	size_t k, const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXXd& theta0,
	size_t max_iter, double tol, bool fix_alpha, bool fix_ab
) : x(x), y(y), z(Eigen::ArrayXd::Ones(x.size())), k(k), theta(theta0), 
	max_iter(max_iter), tol(tol), fix_ab(fix_ab), fix_alpha(fix_alpha), n(x.size())
{
	this->check_init();
	this->gamma = this->get_gamma();
}

Empsvd::EmpsvdCore::EmpsvdCore(
	size_t k, const Eigen::ArrayXd& x, const Eigen::ArrayXd& y,
	size_t max_iter, double tol, bool fix_alpha, bool fix_ab
) : x(x), y(y), z(Eigen::ArrayXd::Ones(x.size())), k(k), 
	max_iter(max_iter), tol(tol), fix_ab(fix_ab), fix_alpha(fix_alpha), n(x.size())
{
	this->theta = this->make_theta0(x, y, k);
	this->check_init();
	this->gamma = this->get_gamma();
}

Empsvd::EmpsvdCore::EmpsvdCore(
	size_t k, const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z, const Eigen::ArrayXXd& theta0,
	size_t max_iter, double tol, bool fix_alpha, bool fix_ab
) : x(x), y(y), z(z), k(k), theta(theta0),
	max_iter(max_iter), tol(tol), fix_ab(fix_ab), fix_alpha(fix_alpha), n(x.size())
{
	this->check_init();
	this->gamma = this->get_gamma();
}

Empsvd::EmpsvdCore::EmpsvdCore(
	size_t k, const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z,
	size_t max_iter, double tol, bool fix_alpha, bool fix_ab
) : x(x), y(y), z(z), k(k),
	max_iter(max_iter), tol(tol), fix_ab(fix_ab), fix_alpha(fix_alpha), n(x.size())
{
	this->theta = this->make_theta0(x, y, z, k);
	this->check_init();
	this->gamma = this->get_gamma();
}

Empsvd::EmpsvdCore::~EmpsvdCore()
{
}

void Empsvd::EmpsvdCore::fit()
{
	double l1 = this->get_loglikelihood();
	double l2;

	for (int i = 0; i < this->max_iter; i++) {
		this->e_step();
		this->m_step();
		this->niter += 1;
		l2 = this->get_loglikelihood();
		if (std::abs(l2 - l1) < this->tol) return;
		l1 = l2;
	}

	throw std::runtime_error("Failed to converge EM algorithm.");
}

void Empsvd::EmpsvdCore::e_step()
{
	this->gamma = this->get_gamma();
}

void Empsvd::EmpsvdCore::m_step()
{
	for (Eigen::Index ik = 0; ik < this->k; ik++) {
		if (!this->fix_ab) {
			this->theta(ik, 2) = this->get_new_bk(ik);
			this->theta(ik, 1) = this->get_new_ak(ik, this->theta(ik, 2));
		}
		this->theta(ik, 3) = this->get_new_sigma2k(ik, this->theta(ik, 1), this->theta(ik, 2));
		if (!this->fix_alpha) {
			this->theta(ik, 4) = this->get_new_alphak(ik);
		}
		this->theta(ik, 5) = this->get_new_lambdak(ik, this->theta(ik, 4));
		this->theta(ik, 0) = this->get_new_omegak(ik);
	}
}

double Empsvd::EmpsvdCore::get_aic()
{
	double l = this->get_loglikelihood();
	return l - (this->k * this->m);
}

double Empsvd::EmpsvdCore::get_bic()
{
	double l = this->get_loglikelihood();
	return l - (0.5 * this->k * this->m * std::log(this->n));
}

double Empsvd::EmpsvdCore::get_loglikelihood()
{
	return this->get_loglikelihood(this->theta);
}

double Empsvd::EmpsvdCore::get_loglikelihood(const Eigen::ArrayXXd& theta)
{
        return (this->get_logsum_pxy(this->get_log_pxy(theta)) * this->z).sum();
}

Eigen::ArrayXXd Empsvd::EmpsvdCore::make_theta0(
	size_t const k, const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXd& z
)
{
	Eigen::ArrayXXd theta0(k, Empsvd::EmpsvdCore::m);

	double const n = z.sum();
	double const y_mean = (y * z).sum() / n;
	double const y_var = ((y - y_mean).square() * z).sum() / n;
	double const x_mean = (x * z).sum() / n;
	double const x_var = ((x - x_mean).square() * z).sum() / n;

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

Eigen::ArrayXXd Empsvd::EmpsvdCore::make_theta0(size_t const k, const Eigen::ArrayXd& x, const Eigen::ArrayXd& y)
{
	return Empsvd::EmpsvdCore::make_theta0(
		x, y, Eigen::ArrayXd::Ones(x.size()), k
	);
}

Eigen::ArrayXXd Empsvd::EmpsvdCore::calc_log_pxy(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXXd& theta)
{
	size_t const n = x.size();
	size_t const k = theta.rows();
	Eigen::ArrayXXd log_pxy(n, k);

	for (Eigen::Index ik = 0; ik < k; ik++) {
		log_pxy.col(ik) = (
			theta(ik, 0) * std::pow(theta(ik, 5), theta(ik, 4)) * x.pow(theta(ik, 4) - 1) /
			(std::sqrt(2 * M_PI * theta(ik, 3)) * std::tgamma(theta(ik, 4)))
			).log() - (
			((y - (theta(ik, 1) * x.pow(theta(ik, 2)))).square() / (2 * theta(ik, 3))) +
				(theta(ik, 5) * x)
				);
	}

	return log_pxy;
}

Eigen::ArrayXXd Empsvd::EmpsvdCore::calc_pxy(const Eigen::ArrayXd& x, const Eigen::ArrayXd& y, const Eigen::ArrayXXd& theta)
{
	return calc_log_pxy(x, y, theta).exp();
}

double Empsvd::EmpsvdCore::digammad(double a)
{
	double const g = 0.57721566490153286061;
	double dig;
	double a2 = std::floor(a);
	double a1 = a - a2;

	if (a1 != 0.) {
		Eigen::ArrayXd c(6);
		c << 0.64493313, -0.20203181, 0.08209433, -0.03591665, 0.01485925, -0.00472050;
		// first, compute the digamma value for 0 < a < 1
		dig = (a1 / (a1 + 1)) - g + (0.5 * std::pow(a1, 7));
		for (int i = 1; i <= c.size(); i++) {
			dig += c(i-1) * (std::pow(a1, i) - std::pow(a1, 7));
		}
		dig -= 1 / a1;
		// then, increment upto a >= 1
		if (a2 > 0.) {
			for (double i = 1.; i <= a2; i++) {
				dig += 1 / (a1 + i - 1);
			}
		}
	}
	else {
		// a1 == 0 means the digamma value is special value
		dig = -g;
		if (a2 > 1.) {
			for (double i = 1; i <= a2 - 1; i++) {
				dig += 1 / i;
			}
		}
	}
	return dig;
}

void Empsvd::EmpsvdCore::check_init()
{
	if (x.size() != n || y.size() != n || z.size() != n) throw std::length_error("Mismatch size between x, y, and z.");
	if (theta.cols() != m) throw std::length_error("Column length of theta must be 6.");
	if (z.sum() <= 0.) throw std::overflow_error("sum(z)=0: This will cause overflow.");
}

Eigen::ArrayXXd Empsvd::EmpsvdCore::get_log_pxy()
{
	return this->get_log_pxy(this->theta);
}

Eigen::ArrayXXd Empsvd::EmpsvdCore::get_log_pxy(const Eigen::ArrayXXd& theta)
{
	return this->calc_log_pxy(this->x, this->y, theta);
}

Eigen::ArrayXXd Empsvd::EmpsvdCore::get_pxy()
{
	return this->get_pxy(this->theta);
}

Eigen::ArrayXXd Empsvd::EmpsvdCore::get_pxy(const Eigen::ArrayXXd& theta)
{
	return this->calc_pxy(this->x, this->y, theta);
}

Eigen::ArrayXd Empsvd::EmpsvdCore::get_logsum_pxy()
{
	return this->get_logsum_pxy(this->get_log_pxy());
}

Eigen::ArrayXd Empsvd::EmpsvdCore::get_logsum_pxy(Eigen::ArrayXXd log_pxy)
{
	Eigen::ArrayXd max_log_pxy = log_pxy.rowwise().maxCoeff();

	for (Eigen::Index ik = 0; ik < this->k; ik++) {
		log_pxy.col(ik) -= max_log_pxy;
	}
	return log_pxy.exp().rowwise().sum().log() + max_log_pxy;
}

Eigen::ArrayXd Empsvd::EmpsvdCore::get_sum_pxy()
{
	return this->get_sum_pxy(this->get_log_pxy());
}

Eigen::ArrayXd Empsvd::EmpsvdCore::get_sum_pxy(Eigen::ArrayXXd log_pxy)
{
	return this->get_logsum_pxy(log_pxy).exp();
}

Eigen::ArrayXXd Empsvd::EmpsvdCore::get_gamma()
{
	return this->get_gamma(this->theta);
}

Eigen::ArrayXXd Empsvd::EmpsvdCore::get_gamma(const Eigen::ArrayXXd& theta)
{
	Eigen::ArrayXXd pxy = this->get_pxy(theta);
	Eigen::ArrayXd sum_pxy = this->get_sum_pxy(pxy.log());

	for (Eigen::Index ik = 0; ik < this->k; ik++) {
		pxy.col(ik) /= sum_pxy;
	}

	return pxy;
}

double Empsvd::EmpsvdCore::get_new_omegak(Eigen::Index ik)
{
	return (this->z * this->gamma.col(ik)).sum() / this->z.sum();
}

double Empsvd::EmpsvdCore::get_new_ak(Eigen::Index ik, double new_bk)
{
	double q1 = (this->z * this->gamma.col(ik) * this->y * this->x.pow(new_bk)).sum();
	double q2 = (this->z * this->gamma.col(ik) * this->x.pow(2 * new_bk)).sum();
	return q1 / q2;
}

double Empsvd::EmpsvdCore::get_new_bk(Eigen::Index ik)
{
	try
	{
		double const offset = 0.001;
		return this->get_new_bk_by_newton(ik, this->theta(ik, 2), this->theta(ik, 2) + offset);
	}
	catch (const std::runtime_error&)
	{
		size_t const max_retry = 30;
		Eigen::ArrayXd bk1 = (Eigen::ArrayXd::Random(max_retry) + 1.0) * 0.5; // [0:1] random number
		for (Eigen::Index i = 0; i < max_retry; i++) {
			try
			{
				return this->get_new_bk_by_newton(ik, this->theta(ik, 2), bk1(i));
			}
			catch (const std::runtime_error& e)
			{
				continue;
			}
		}
		throw std::runtime_error("Failed to find an optimal new bk.");
	}
}

double Empsvd::EmpsvdCore::get_new_sigma2k(Eigen::Index ik, double new_ak, double new_bk)
{
	double q1 = (this->z * this->gamma.col(ik) * (this->y - (new_ak * this->x.pow(new_bk))).square()).sum();
	double q2 = (this->z * this->gamma.col(ik)).sum();
	return q1 / q2;
}

double Empsvd::EmpsvdCore::get_new_alphak(Eigen::Index ik)
{
	return this->get_new_alphak_by_invdigamma(ik, this->theta(ik, 4));
}

double Empsvd::EmpsvdCore::get_new_lambdak(Eigen::Index ik, double new_alphak)
{
	double q1 = (this->z * this->gamma.col(ik)).sum();
	double q2 = (this->z * this->gamma.col(ik) * x).sum();
	return new_alphak * q1 / q2;
}

double Empsvd::EmpsvdCore::bkdot(Eigen::Index ik, double bk)
{
	double ak = this->get_new_ak(ik, bk);
	return (
		this->z * this->gamma.col(ik) * this->x.log() *
		((this->y * this->x.pow(bk)) - (ak * this->x.pow(2 * bk)))
		).sum();
}

double Empsvd::EmpsvdCore::get_new_bk_by_newton(Eigen::Index ik, double bk0, double bk1, size_t max_iter, double tol)
{
	double fbk0 = this->bkdot(ik, bk0);
	double fbk1 = this->bkdot(ik, bk1);
	double fbk2, bk2;

	for (int i = 0; i < max_iter; i++) {
		bk2 = ((bk0 * fbk1) - (bk1 * fbk0)) / (fbk1 - fbk0);
		fbk2 = this->bkdot(ik, bk2);

		if (std::abs(fbk2) < tol) return bk2;

		bk0 = bk1;
		fbk0 = fbk1;
		bk1 = bk2;
		fbk1 = fbk2;
	}
	throw std::runtime_error("Failed to find optimal bk with newton method.");
}

double Empsvd::EmpsvdCore::get_new_alphak_by_invdigamma(Eigen::Index ik, double alphak0, size_t max_iter, double tol)
{
	double q1 = (this->z * this->gamma.col(ik)).sum();
	double logx_mean = (this->z * this->gamma.col(ik) * this->x.log()).sum() / q1;
	double x_mean = (this->z * this->gamma.col(ik) * this->x).sum() / q1;
	double y = std::log(x_mean) - logx_mean;

	double dig, alphak1;
	
	for (int i = 0; i < max_iter; i++) {
		dig = this->digammad(alphak0);
		alphak1 = alphak0 * (std::log(alphak0) - dig) / y;
		if (std::abs(alphak1 - alphak0) < tol) return alphak1;
		alphak0 = alphak1;
	}
	throw std::runtime_error("Failed to find an optimal alpha.");
}

