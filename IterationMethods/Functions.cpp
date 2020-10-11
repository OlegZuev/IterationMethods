#include "Functions.h"
#include <utility>
#include "Vector.h"

double find_transition_matrix_norm(const Vector& prev_x, const Vector& x, const Vector& next_x) {
	return (next_x - x).find_third_norm() / (x - prev_x).find_third_norm();
}

double error_estimate(double q, const Vector& x, const Vector& next_x) {
	return (next_x - x).find_third_norm() * q / (1 - q);
}

void print_header(std::ostream& ostr){
	ostr << "                        | Discrepancy |  Error    |" << std::endl;
	ostr << "Itr |    tau   |  q     | Norm        |  Estimate |     x[1]      x[2]      x[3]      x[4]" << std::endl;
}