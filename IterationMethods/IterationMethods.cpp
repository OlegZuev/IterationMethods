#include <iostream>
#include <string>
#include <fstream>
#include "Matrix.h"
using namespace std;

int main() {
	string variant = "input16b";
	ifstream fin("../" + variant + ".txt");

	int n;
	fin >> n;
	Matrix matr(n);
	matr.input(fin);
	Vector b(n);
	b.input(fin);

	matr.print(cout);

	cout << matr.get_third_norm() << std::endl;

	fin.close();

	ofstream fout("../" + variant + "_output" + ".txt");
	ostream& out = cout;	

	Vector result = matr.conjugate_gradient_method(b, out);
	result.print(std::cout);

	fout.close();
}
