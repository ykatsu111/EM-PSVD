#include "EMPSVD.h"
#include <Eigen/Dense>

EMPSVD::EMPSVD(Eigen::ArrayXd& x, Eigen::ArrayXd& y, int k, Eigen::ArrayXXd theta0)
{
	this->x = x;
	this->y = y;
	this->k = k;
	this->theta = theta0;

	this->n = x.size();
}

EMPSVD::EMPSVD(Eigen::ArrayXd& x, Eigen::ArrayXd& y, int k)
{
	Eigen::ArrayXXd theta = this->make_theta0();
	EMPSVD(x, y, k, theta);
}

EMPSVD::~EMPSVD()
{
}

Eigen::ArrayXXd EMPSVD::make_theta0(Eigen::ArrayXd& const x, Eigen::ArrayXd& const y, int const k, int const m)
{
	Eigen::ArrayXXd theta0(k, m);
	Eigen::ArrayXd const ones = Eigen::ArrayXd::Ones(m);

	int const n = x.size();
	double const y_mean = y.mean();
	double const y_var = (y - y_mean).square().sum() / n;
	double const x_mean = x.mean();
	double const x_var = (x - x_mean).square().sum() / n;

	theta0.col(0) = 1. / k;
	theta0.col(3) = y_var;
	theta0.col(4) = std::pow(x_mean, 2) / x_var;
	theta0.col(5) = x_mean / x_var;

	double a_upp = y.maxCoeff() / x_mean;
	double a_low = y_mean / 4;

	if (k > 1) {
		Eigen::ArrayXd lins = Eigen::ArrayXd::LinSpaced(k, 0., 1.);
		theta0.col(1) = (lins * std::log(a_upp / a_low) + std::log(a_low)).exp();
		theta0.col(1) = lins;
	}
	else {
		theta0.col(1) = std::exp(std::log(a_upp * a_low) / 2.);
		theta0.col(1) = 0.5;
	}

	return theta0;
}
