#include "EmpsvdBinCore.h"

Empsvd::EmpsvdBinCore::EmpsvdBinCore(
	Eigen::ArrayXXd x, Eigen::ArrayXXd y, Eigen::ArrayXXd z, size_t k, Eigen::ArrayXXd theta0, 
	size_t max_iter, double tol, bool fix_alpha, bool fix_ab
) : Empsvd::EmpsvdCore(
	Eigen::Map<Eigen::ArrayXd>(x.data(), x.size()), Eigen::Map<Eigen::ArrayXd>(y.data, y.size()), k, theta0,
	max_iter, tol, fix_alpha, fix_ab
), z(Eigen::Map<Eigen::ArrayXd>(z.data(), z.size()))
{
}

Empsvd::EmpsvdBinCore::EmpsvdBinCore(
	Eigen::ArrayXXd x, Eigen::ArrayXXd y, Eigen::ArrayXXd z, size_t k, 
	size_t max_iter, double tol, bool fix_alpha, bool fix_ab
) : Empsvd::EmpsvdCore(
	Eigen::Map<Eigen::ArrayXd>(x.data(), x.size()), Eigen::Map<Eigen::ArrayXd>(y.data, y.size()), k, 
	this->make_theta0(x, y, z, k, this->m),
	max_iter, tol, fix_alpha, fix_ab
), z(Eigen::Map<Eigen::ArrayXd>(z.data(), z.size()))
{
}
