#include <iostream>
#include <string>
#include <fstream>
#include "Matrix.h"
using namespace std;

const bool PRIN_TO_FILE = false;

int main() {
	string variant = "input1";
	ifstream fin("../" + variant + ".txt");

	int n;
	fin >> n;
	Matrix matr(n);
	matr.input(fin);
	Vector b(n);
	b.input(fin);
	fin.close();

	ofstream fout("../" + variant + "_output" + ".txt");
	ostream& out = PRIN_TO_FILE ? fout: cout;

	out << "Variand " << variant;
	out << "b" << endl;
	b.print(out);

	out << "A" << endl;
	matr.print(cout);

	cout << "Third norm of a matrix " << matr.get_third_norm() << std::endl << std::endl;

	Vector result(n);
	result = matr.simple_iteration_method(b, out);
	result.print(std::cout);
	out << endl;
	result = matr.gradient_steepest_descent_method(b, out);
	result.print(std::cout);
	out << endl;
	result = matr.sor_method(b, out);
	result.print(std::cout);
	out << endl;
	result = matr.conjugate_gradient_method(b, out);
	result.print(std::cout);
	out << endl;

	fout.close();
}
