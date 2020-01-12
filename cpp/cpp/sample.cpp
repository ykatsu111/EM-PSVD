#include "EmpsvdCore.h"
#include "csv.h"
#include <iostream>
#include <Eigen/Dense>

constexpr int COLUMNS = 2;
constexpr int ROWS = 900;

int main() {
	io::CSVReader<COLUMNS> in("../sample_data/psvd.csv");
	Eigen::ArrayXd x(ROWS), y(ROWS);
	double xi, yi;
	Eigen::Index i = 0;
	
	while (in.read_row(xi, yi)) {
		x(i) = xi;
		y(i) = yi;
		i++;
	}

	EmpsvdCore em(x, y, 2);
	std::cout << em.theta << std::endl;
	em.fit();
	std::cout << em.theta << std::endl;

	return 0;
}