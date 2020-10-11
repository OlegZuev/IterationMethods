#include <iostream>
#include <string>
#include <fstream>
#include "Matrix.h"
#include <iomanip>
using namespace std;

const bool PRIN_TO_FILE = true;

int main() {
	string variant = "input7b";

	int n;

	// Input n, A, b
	ifstream fin("../" + variant + ".txt");
	fin >> n;
	Matrix matr(n);
	Vector b(n);
	matr.input(fin);
	b.input(fin);
	fin.close();

	ofstream fout("../" + variant + "_output" + ".txt");
	ostream& out = PRIN_TO_FILE ? fout : cout;

	out << "Variant " << variant << endl;
	out << "b" << endl;
	b.print(out);

	out << "A" << endl;
	matr.print(out);

	out << "Third norm of a matrix " << matr.get_third_norm() << std::endl << std::endl;

	Vector result(n);

	result = matr.simple_iteration_method(b, out);
	result.print(out);
	out << endl;

	result = matr.gradient_steepest_descent_method(b, out);
	result.print(out);
	out << endl;

	result = matr.sor_method(b, out);
	result.print(out);
	out << endl;

	result = matr.conjugate_gradient_method(b, out);
	result.print(out);
	out << endl;

	double condition_number = matr.get_third_norm() * matr.get_inverse_matrix().get_third_norm();
	out << "Condition number: " << condition_number << endl;
	out << "Theoretical estimation of the number of iterations: " << endl;
	out << setprecision(0) << "Simple iteration method: " << round(log(1 / EPS) / 2 * condition_number) << endl;
	out << setprecision(0) << "Gradient Steepest descent method: " << round(log(1 / EPS) / 2 * condition_number) << endl;
	out << setprecision(0) << "SOR method: " << round(log(1 / EPS) / 4 * sqrt(condition_number)) << endl;
	out << setprecision(0) << "Conjugate gradient method: " << round(log(2 / EPS) / 2 * sqrt(condition_number)) << endl;

	fout.close();
}
