#pragma once
#include <iostream>

class Vector {
	int size{};
	double* arr{};

public:
	static constexpr double eps = 1E-4;

	Vector() = delete;
	explicit Vector(int n);
	Vector(const Vector& other);
	Vector& operator=(Vector other);
	Vector(Vector&& other) noexcept;
	~Vector();

	friend void swap(Vector& lhs, Vector& rhs) noexcept;

	double& operator[](int index) const;
	void input(std::istream& fin) const;
	void print(std::ostream& ostr) const;
	double find_third_norm() const;
	Vector operator-(const Vector& other) const;
	Vector operator+(const Vector& other) const;
	Vector operator*(const double number) const;
};