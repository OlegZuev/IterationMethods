#pragma once
#include <iostream>
#include "Vector.h"
#define PRECISION 7
#define EPS 1E-4
#define EPS_COMPUTING_W 1E-2

class Matrix {
	int size{};
	double** arr{};

public:
	Matrix() = delete;
	explicit Matrix(int n);
	Matrix(const Matrix& other);
	Matrix& operator=(Matrix other);
	Matrix(Matrix&& other) noexcept;
	~Matrix();

	friend void swap(Matrix& lhs, Matrix& rhs) noexcept;

	void input(std::istream& fin) const;
	void print(std::ostream& ostr) const;
	double*& operator[](int index) const;
	Matrix operator*(const Matrix& other) const;
	Vector operator*(const Vector& vec) const;
	double get_first_norm() const;
	double get_second_norm() const;
	double get_third_norm() const;
	Matrix get_diagonal_matrix() const;
	double max_not_diagonal(int& i, int& j) const;
	Vector simple_iteration_method(const Vector& b, std::ostream& ostr) const;
	Vector gradient_steepest_descent_method(const Vector& b, std::ostream& ostr) const;
	void find_next_x(double w, const Vector& x, Vector& next_x, const Vector& b) const;
	double find_optimal_w(const Vector& b, std::ostream& ostr) const;
	Vector sor_method(const Vector& b, std::ostream& ostr) const;
	Vector conjugate_gradient_method(const Vector& b, std::ostream& ostr) const;
	void find_LU(Matrix& matrix_PA, Matrix& matrix_U, Matrix& matrix_L, Vector& indexes, int& permutation_count) const;
	void find_inverse_matrix(Matrix& matrix_inverseA, const Matrix& matrix_U, const Matrix& matrix_L, const Vector& indexes) const;
	Matrix get_inverse_matrix() const;
	Vector get_x_with_LU(const Vector& b) const;
	Vector find_x(const Vector& b, const Matrix& matrix_U, const Matrix& matrix_L, const Vector& indexes) const;
};