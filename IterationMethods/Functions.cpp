#include "Functions.h"
#include <utility>
#include "Vector.h"

double find_third_norm_vector(double* x, int n) {
	double res = 0;
	for (int i = 0; i < n; i++) {
		res += x[i] * x[i];
	}

	return sqrt(res);
}

double find_transition_matrix_norm(const int itr, const Vector& prev_x, const Vector& x, const Vector& next_x) {
	return (next_x - x).find_third_norm() / (x - prev_x).find_third_norm();
}

double error_estimate(double q, const Vector& prev_x, const Vector& x) {
	return (prev_x - x).find_third_norm() * q / (1 - q);
}