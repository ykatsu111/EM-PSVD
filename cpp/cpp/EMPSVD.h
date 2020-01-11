#pragma once

#include<Eigen/Dense>

class EMPSVD
{
public:
	EMPSVD(
		Eigen::ArrayXd& x,
		Eigen::ArrayXd& y,
		int k,
		Eigen::ArrayXXd theta0
	);

	EMPSVD(
		Eigen::ArrayXd& x,
		Eigen::ArrayXd& y,
		int k
	);
	void fit();
	double get_aic();
	double get_bic();
	double get_loglikelihood();

	~EMPSVD();

	Eigen::ArrayXd x, y;
	int k, n;
	int const m = 6;
	Eigen::ArrayXXd theta;

private:
	static Eigen::ArrayXXd make_theta0(
		Eigen::ArrayXd& const x,
		Eigen::ArrayXd& const y,
		int const k,
		int const m
	);


};

