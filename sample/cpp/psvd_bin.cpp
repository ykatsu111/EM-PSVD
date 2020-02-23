#include "EmpsvdBinCore.h"
#include "csv.h"
#include <iostream>
#include <Eigen/Dense>

constexpr int COLUMNS = 3;
constexpr int ROWS = 100 * 60;

int main() {
	io::CSVReader<COLUMNS> in("psvd_bin.csv");
	Eigen::ArrayXd x(ROWS), y(ROWS), z(ROWS);
	double xi, yi, zi;
	Eigen::Index i = 0;
	
	while (in.read_row(xi, yi, zi)) {
		x(i) = xi;
		y(i) = yi;
		z(i) = zi;
		i++;
	}

	Empsvd::EmpsvdBinCore em(x, y, z, 2, 1000, 1e-5);

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
