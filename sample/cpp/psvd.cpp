#include "EmpsvdCore.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <Eigen/Dense>

constexpr int COLUMNS = 2;
constexpr int ROWS = 900;

int main() {
	std::ifstream ifs("../data/psvd.csv");
	Eigen::ArrayXd x(ROWS), y(ROWS);
	Eigen::Index i = 0;
	std::string buf, field;
	
	if (!ifs.is_open()) {
		std::cerr << "No file is opened." << std::endl;
		return 0;
	}
	std::getline(ifs, buf);  // skip header
	while (std::getline(ifs, buf)) {
		std::stringstream ss(buf);
		std::getline(ss, field, ','); x(i) = std::stod(field);
		std::getline(ss, field, ','); y(i) = std::stod(field);
		i++;		
	}

	Empsvd::EmpsvdCore em(2, x, y, 1000, 1e-5);

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
