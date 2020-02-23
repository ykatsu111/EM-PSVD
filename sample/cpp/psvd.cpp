#include "EmpsvdCore.h"
#include "csv.h"
#include <iostream>
#include <Eigen/Dense>

constexpr int COLUMNS = 2;
constexpr int ROWS = 900;

int main() {
	io::CSVReader<COLUMNS> in("psvd.csv");
	Eigen::ArrayXd x(ROWS), y(ROWS);
	double xi, yi;
	Eigen::Index i = 0;
	
	while (in.read_row(xi, yi)) {
		x(i) = xi;
		y(i) = yi;
		i++;
	}

	Empsvd::EmpsvdCore em(x, y, 2, 1000, 1e-5);

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
