/*
EM-PSVD
(C) Y. Katsuyama (2020)
See https://humet.sci.hokudai.ac.jp/~meteo/product-e.html
*/

#include "EmpsvdCore.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <Eigen/Dense>

constexpr int ROWS = 100 * 60;

int main() {
	std::ifstream ifs("../data/psvd_bin.csv");
	Eigen::ArrayXd x_in(ROWS), y_in(ROWS), z_in(ROWS);
	Eigen::Index i = 0;
	std::string buf, field;
	
	if (!ifs.is_open()) {
		std::cerr << "No file is opened." << std::endl;
		return 0;
	}
	std::getline(ifs, buf);  // skip header
	while (std::getline(ifs, buf)) {
		std::stringstream ss(buf);
		std::getline(ss, field, ','); x_in(i) = std::stod(field);
		std::getline(ss, field, ','); y_in(i) = std::stod(field);
		std::getline(ss, field, ','); z_in(i) = std::stod(field);
		i++;
	}

	// extract elements where z > 0
	Eigen::ArrayXd exist_data = (z_in.array() > 0.).select(					     
		Eigen::ArrayXd::Constant(z_in.size(), 1),
		Eigen::ArrayXd::Constant(z_in.size(), 0)
	);
	size_t n = exist_data.count();
	Eigen::ArrayXd x(n), y(n), z(n);
	Eigen::Index j = 0;
	for (i = 0; i < x_in.size(); i++) {
	        if (exist_data(i)) {
		        x(j) = x_in(i);
			y(j) = y_in(i);
			z(j) = z_in(i);
			j++;
	        }
	}

	Empsvd::EmpsvdCore em(2, x, y, z, 1000, 1e-5);

	std::cout << "---Initial States---" << std::endl;
	std::cout << "omega   a     b     sigma2     alpha(mu+1)    lambda" << std::endl;
	std::cout << em.theta << std::endl;
	std::cout << "Loglikelihood: " << em.get_loglikelihood() << std::endl;
	std::cout << "AIC: " << em.get_aic() << std::endl;
	std::cout << "BIC: " << em.get_bic() << std::endl;

	em.fit();
	std::cout << "---Fitting Result----" << std::endl;
	std::cout << "omega   a     b     sigma2     alpha(mu+1)    lambda" << std::endl;
	std::cout << em.theta << std::endl;
	std::cout << "Loglikelihood: " << em.get_loglikelihood() << std::endl;
	std::cout << "AIC: " << em.get_aic() << std::endl;
	std::cout << "BIC: " << em.get_bic() << std::endl;
	std::cout << "Total iteration: " << em.niter << std::endl;

	return 0;
}
