#include "Matrix.h"
#include "Functions.h"
#include <utility>
#include <iomanip>

/**
 * Constructor with allocating memory for matrix
 */
Matrix::Matrix(int n) {
	size = n;
	arr = new double* [size];
	for (int i = 0; i < size; ++i) {
		arr[i] = new double[size] {0};
	}
}

/**
 * Copy constructor
 */
Matrix::Matrix(const Matrix& other) : Matrix(other.size) {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			arr[i][j] = other.arr[i][j];
		}
	}
}

/**
 * Copy assignment operator
 */
Matrix& Matrix::operator=(Matrix other) {
	swap(*this, other);

	return *this;
}

/**
 * Move constructor
 */
Matrix::Matrix(Matrix&& other) noexcept {
	swap(*this, other);
}

/**
 * Destructor
 */
Matrix::~Matrix() {
	for (int i = 0; i < size; ++i) {
		delete[] arr[i];
	}

	delete[] arr;
}

/**
 * Swap
 */
void swap(Matrix& lhs, Matrix& rhs) noexcept {
	using std::swap;

	swap(lhs.size, rhs.size);
	swap(lhs.arr, rhs.arr);
}

/**
 * Input matrix from "fin" stream
 */
void Matrix::input(std::istream& fin) const {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			fin >> arr[i][j];
		}
	}
}

/**
 * Print matrix into "ostr" stream
 */
void Matrix::print(std::ostream& ostr) const {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			ostr << arr[i][j] << " ";
		}

		ostr << std::endl;
	}
}

/**
 * Indexer
 */
double*& Matrix::operator[](const int index) const {
	return arr[index];
}

/**
 * Operator for multiply matrix: like A*B
 */
Matrix Matrix::operator*(const Matrix& other) const {
	Matrix new_matrix(size);
	for (int i = 0; i < size; ++i) {
		for (int k = 0; k < size; ++k) {
			for (int j = 0; j < size; ++j) {
				new_matrix[i][k] += (*this)[i][j] * other[j][k];
			}
		}
	}

	return new_matrix;
}

/**
 * Operator for multiply matrix to vector: like A*x
 */
Vector Matrix::operator*(const Vector& vec) const {
	Vector result(size);
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			result[i] += (*this)[i][j] * vec[j];
		}
	}

	return result;
}

/**
 * First norm of matrix
 */
double Matrix::get_first_norm() const {
	double result = -1;
	for (int i = 0; i < size; ++i) {
		double sum = 0;
		for (int j = 0; j < size; ++j) {
			sum += fabs((*this)[i][j]);
		}

		if (sum > result) {
			result = sum;
		}
	}

	return result;
}

/**
 * Second norm of matrix
 */
double Matrix::get_second_norm() const {
	double result = -1;
	for (int i = 0; i < size; ++i) {
		double sum = 0;
		for (int j = 0; j < size; ++j) {
			sum += fabs((*this)[j][i]);
		}

		if (sum > result) {
			result = sum;
		}
	}

	return result;
}

/**
 * Third norm of matrix
 */
double Matrix::get_third_norm() const {
	Matrix tmp_matrix(size);
	// A * A^T
	for (int i = 0; i < size; ++i) {
		for (int k = 0; k < size; ++k) {
			for (int j = 0; j < size; ++j) {
				tmp_matrix[i][k] += (*this)[j][i] * (*this)[j][k];
			}
		}
	}

	// Diagonal matrix for A * A^T
	Matrix diagonal_matrix = tmp_matrix.get_diagonal_matrix();
	// Finding max diagonal element
	double result = fabs(diagonal_matrix[0][0]);
	for (int i = 1; i < size; ++i) {
		if (result < fabs(diagonal_matrix[i][i])) {
			result = fabs(diagonal_matrix[i][i]);
		}
	}

	return sqrt(result);
}

/**
 * Get diagonal matrix of eigenvalues
 */
Matrix Matrix::get_diagonal_matrix() const {
	int im = 0;
	int jm = 0;
	// Copying matrix
	Matrix new_matrix = (*this);

	Matrix tmp_matrix(size);

	while (new_matrix.max_not_diagonal(im, jm) > EPS) {
		double alpha = atan(2 * new_matrix[im][jm] / (new_matrix[im][im] - new_matrix[jm][jm])) / 2;

		Matrix matrix_T(size);
		for (int i = 0; i < matrix_T.size; ++i) {
			matrix_T[i][i] = 1;
		}

		matrix_T[im][im] = cos(alpha);
		matrix_T[jm][jm] = cos(alpha);
		matrix_T[im][jm] = -sin(alpha);
		matrix_T[jm][im] = sin(alpha);

		// T^T * A
		for (int i = 0; i < size; ++i) {
			for (int k = 0; k < size; ++k) {
				tmp_matrix[i][k] = 0;
				for (int j = 0; j < size; ++j) {
					tmp_matrix[i][k] += matrix_T[j][i] * new_matrix[j][k];
				}
			}
		}

		new_matrix = tmp_matrix * matrix_T;
	}

	return new_matrix;
}

/**
 * Find max element in matrix. This element must not diagonal
 */
double Matrix::max_not_diagonal(int& i, int& j) const {
	double res = -1;
	for (int k = 0; k < size; ++k) {
		for (int p = 0; p < size; ++p) {
			if (k == p)
				continue;

			if (res < fabs((*this)[k][p])) {
				res = fabs((*this)[k][p]);
				i = k;
				j = p;
			}
		}
	}

	return res;
}

/**
 * Simple iteration method for solving system of linear equations
 */
Vector Matrix::simple_iteration_method(const Vector& b, std::ostream& ostr) const {
	ostr << "Simple iteration method" << std::endl;
	print_header(ostr);
	Vector x = b;
	Vector prev_x(size);
	Vector next_x(size);
	int itr = 1;
	const double tau = 1.8 / get_third_norm();
	double discrepancy_norm;
	do {
		// find next x
		for (int i = 0; i < size; ++i) {
			double sum = 0;
			for (int j = 0; j < size; ++j) {
				sum += (*this)[i][j] * x[j];
			}

			next_x[i] = x[i] + tau * (b[i] - sum);
		}

		// ||b - A * x(k+1)||
		discrepancy_norm = (b - (*this) * next_x).find_third_norm();
		double q = find_transition_matrix_norm(prev_x, x, next_x);
		double e = error_estimate(q, x, next_x);

		prev_x = x;
		x = next_x;
		ostr << std::fixed << std::setprecision(PRECISION) << itr << " | " << tau << " | " << q << " | " <<
			discrepancy_norm << " | " << e << " | " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;
		itr++;
	} while (discrepancy_norm > EPS);

	return x;
}

/**
 * Gradient steepest descent method for solving system of linear equations
 */
Vector Matrix::gradient_steepest_descent_method(const Vector& b, std::ostream& ostr) const {
	ostr << "Gradient steepest descent method" << std::endl;
	print_header(ostr);
	Vector x = b;
	Vector prev_x(size);
	Vector next_x(size);
	int itr = 1;
	double tau;
	double discrepancy_norm;

	do {
		// find new tau
		Vector r = (*this) * next_x - b;
		tau = (r * r) / (((*this) * r) * r);

		// find next x
		for (int i = 0; i < size; ++i) {
			double sum = 0;
			for (int j = 0; j < size; ++j) {
				sum += (*this)[i][j] * x[j];
			}

			next_x[i] = x[i] + tau * (b[i] - sum);
		}

		// ||b - A * x(k+1)||
		discrepancy_norm = (b - (*this) * next_x).find_third_norm();
		double q = find_transition_matrix_norm(prev_x, x, next_x);
		double e = error_estimate(q, x, next_x);

		prev_x = x;
		x = next_x;
		ostr << std::fixed << std::setprecision(PRECISION) << itr << " | " << tau << " | " << q << " | " <<
			discrepancy_norm << " | " << e << " | " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;
		itr++;
	} while (discrepancy_norm > EPS);

	return x;
}

/**
 * finding next vector x for SOR
 */
void Matrix::find_next_x(double w, const Vector& x, Vector& next_x, const Vector& b) const {
	for (int i = 0; i < size; i++) {
		double sum1 = 0;
		double sum2 = 0;

		for (int j = 0; j <= i - 1; j++) {
			sum1 += (*this)[i][j] * next_x[j];
		}

		for (int j = i + 1; j < size; j++) {
			sum2 += (*this)[i][j] * x[j];
		}

		double x_zeydel = (1 / (*this)[i][i]) * (b[i] - sum1 - sum2);

		next_x[i] = x[i] + w * (x_zeydel - x[i]);
	}
}

/**
 * finding optimal omega (w) for SOR method
 */
double Matrix::find_optimal_w(const Vector& b, std::ostream& ostr) const {
	double w = 0.1;
	double optimal_w = 0.1;
	int min_itr = INT_MAX;
	double discrepancy_norm;

	ostr << "SOR method - search optimal w" << std::endl;
	while (w < 2) {
		Vector next_x(size);
		Vector x = b;
		int itr = 1;
		do {
			find_next_x(w, x, next_x, b);

			// ||b - A * x(k+1)||
			discrepancy_norm = (b - (*this) * next_x).find_third_norm();

			x = next_x;
			itr++;
		} while (discrepancy_norm > EPS_COMPUTING_W);

		if (itr < min_itr) {
			min_itr = itr;
			optimal_w = w;
		}

		ostr << "w = " << w << " Itr = " << itr << std::endl;
		w += 0.1;
	}

	ostr << std::endl;
	ostr << "w = " << optimal_w << " ItrMin = " << min_itr << std::endl;
	return optimal_w;
}

/**
 * Successive over-relaxation method for solving system of linear equations
 */
Vector Matrix::sor_method(const Vector& b, std::ostream& ostr) const {
	int itr = 1;
	double w = find_optimal_w(b, ostr);
	double discrepancy_norm;
	Vector x = b;
	Vector prev_x(size);
	Vector next_x(size);

	ostr << "SOR method" << std::endl;
	print_header(ostr);

	do {
		find_next_x(w, x, next_x, b);

		// ||b - A * x(k+1)||
		discrepancy_norm = (b - (*this) * next_x).find_third_norm();
		double q = find_transition_matrix_norm(prev_x, x, next_x);
		double e = error_estimate(q, x, next_x);

		prev_x = x;
		x = next_x;
		ostr << std::fixed << std::setprecision(PRECISION) << itr << " | " << w << " | " << q << " | " <<
			discrepancy_norm << " | " << e << " | " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;
		itr++;
	} while (discrepancy_norm > EPS);

	return x;
}

/**
 * Conjugate gradient method for solving system of linear equations
 */
Vector Matrix::conjugate_gradient_method(const Vector& b, std::ostream& ostr) const {
	int itr = 1;
	double discrepancy_norm;
	Vector x = b; // x(k)
	Vector prev_x(size); // x(k - 1)
	Vector next_x(size); // x(k + 1)
	Vector prev_r(size); // 
	double prev_rr = 0;
	double tau; // tau(k + 1)
	double prev_tau = 0; // tau
	double alpha = 1;

	ostr << "Conjugate gradient method" << std::endl;
	print_header(ostr);

	do {
		// computing Ax - b
		Vector r = (*this) * x - b;
		double rr = r * r;

		// find new tau
		tau = (r * r) / (((*this) * r) * r);

		// For first iteration alpha = 1
		if (itr != 1) {
			alpha = 1 / (1 - tau * rr / prev_tau / alpha / prev_rr);
		}

		// Computing new x
		for (int i = 0; i < size; i++) {
			next_x[i] = alpha * x[i] + (1 - alpha) * prev_x[i] - tau * alpha * r[i];
		}

		discrepancy_norm = (b - (*this) * next_x).find_third_norm();
		double q = find_transition_matrix_norm(prev_x, x, next_x);
		double e = error_estimate(q, x, next_x);

		prev_x = x;
		x = next_x;
		prev_tau = tau;
		prev_r = r;
		prev_rr = rr;
		ostr << std::fixed << std::setprecision(PRECISION) << itr << " | " << tau << " | " << q << " | " <<
			discrepancy_norm << " | " << e << " | " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " | " << alpha << std::endl;
		itr++;
	} while (discrepancy_norm > EPS);

	return x;
}

/**
 * Find LU for matrix with permutations in indexes
 */
void Matrix::find_LU(Matrix& matrix_PA, Matrix& matrix_U, Matrix& matrix_L, Vector& indexes, int& permutation_count) const {
	matrix_PA = (*this);
	matrix_U = (*this);
	permutation_count = 0;
	int rank = size;
	for (int i = 0; i < size; ++i) {
		double column_max = matrix_U[i][i];
		int k = i;
		for (int j = k + 1; j < size; ++j) {
			if (fabs(column_max) < fabs(matrix_U[j][i])) {
				column_max = matrix_U[j][i];
				k = j;
			}
		}

		if (fabs(matrix_PA[k][i]) < EPS) {
			rank--;
		}

		if (i != k) {
			permutation_count++;
		}

		std::swap(matrix_U[i], matrix_U[k]);
		std::swap(indexes[i], indexes[k]);
		std::swap(matrix_PA[i], matrix_PA[k]);

		double multiplier = matrix_U[i][i];
		for (int j = 0; j < size; ++j) {
			matrix_U[i][j] = matrix_U[i][j] / multiplier;
		}

		for (int b = i + 1; b < size; ++b) {
			multiplier = matrix_U[b][i];
			for (int j = i; j < size; ++j) {
				matrix_U[b][j] = matrix_U[b][j] - multiplier * matrix_U[i][j];
			}
		}

		for (int j = 0; j <= i; ++j) {
			double vector_mult_result = 0;
			for (int k2 = 0; k2 < i; ++k2) {
				vector_mult_result += matrix_L[i][k2] * matrix_U[k2][j];
			}

			matrix_L[i][j] = matrix_PA[i][j] - vector_mult_result;
		}
	}
}

/**
 * Find inverse matrix
 */
void Matrix::find_inverse_matrix(Matrix& matrix_inverseA, const Matrix& matrix_U, const Matrix& matrix_L, const Vector& indexes) const {
	Vector y(size);

	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			double sum = 0;
			for (int k = 0; k < j; ++k) {
				sum += matrix_L[j][k] * y[k];
			}

			y[j] = ((i == indexes[j] ? 1 : 0) - sum) / matrix_L[j][j];
		}

		for (int j = size - 1; j >= 0; --j) {
			double sum = 0;
			for (int k = j + 1; k < size; ++k) {
				sum += matrix_U[j][k] * matrix_inverseA[k][i];
			}

			matrix_inverseA[j][i] = y[j] - sum;
		}
	}
}

/**
 * Get inverse matrix (via Lu decomposition)
 */
Matrix Matrix::get_inverse_matrix() const {
	int permutation_count;
	Matrix matrix_PA(size);
	Matrix matrix_L(size);
	Matrix matrix_U(size);
	Matrix matrix_inverseA(size);
	Vector indexes(size);
	for (int i = 0; i < size; ++i) {
		indexes[i] = i;
	}

	find_LU(matrix_PA, matrix_U, matrix_L, indexes, permutation_count);
	find_inverse_matrix(matrix_inverseA, matrix_U, matrix_L, indexes);
	return matrix_inverseA;
}

/**
 * Solving linear system (via LU)
 */
Vector Matrix::get_x_with_LU(const Vector& b) const {
	int permutation_count;
	Matrix matrix_PA(size);
	Matrix matrix_L(size);
	Matrix matrix_U(size);
	Matrix matrix_inverseA(size);
	Vector indexes(size);
	for (int i = 0; i < size; ++i) {
		indexes[i] = i;
	}

	find_LU(matrix_PA, matrix_U, matrix_L, indexes, permutation_count);
	return find_x(b, matrix_U, matrix_L, indexes);
}

/**
 * Solving linear system with L and U
 */
Vector Matrix::find_x(const Vector& b, const Matrix& matrix_U, const Matrix& matrix_L, const Vector& indexes) const {
	Vector y(size);
	Vector x = b;
	for (int j = 0; j < size; ++j) {
		double sum = 0;
		for (int k = 0; k < j; ++k) {
			sum += matrix_L[j][k] * y[k];
		}

		y[j] = (b[indexes[j]] - sum) / matrix_L[j][j];
	}

	for (int j = size - 1; j >= 0; --j) {
		double sum = 0;
		for (int k = j + 1; k < size; ++k) {
			sum += matrix_U[j][k] * x[k];
		}

		x[j] = y[j] - sum;
	}

	return x;
}
