#include "EmpsvdCore.h"
#include "csv.h"
#include <iostream>
#include <Eigen/Dense>

constexpr int COLUMNS = 3;
constexpr int ROWS = 100 * 60;

int main() {
	io::CSVReader<COLUMNS> in("psvd_bin.csv");
	Eigen::ArrayXd x_in(ROWS), y_in(ROWS), z_in(ROWS);
	double xi, yi, zi;
	Eigen::Index i = 0;
	
	while (in.read_row(xi, yi, zi)) {
		x_in(i) = xi;
		y_in(i) = yi;
		z_in(i) = zi;
		i++;
	}
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

	Empsvd::EmpsvdCore em(x, y, z, 2, 1000, 1e-5);

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
