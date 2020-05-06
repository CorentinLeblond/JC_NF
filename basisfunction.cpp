#include "basisfunction.hpp"
#include<cmath>

basis_functions::basis_functions(int order) :
base_order(order)
{};

Poly_Laguerre::Poly_Laguerre(int order) :
	basis_functions(order)
{};

double Poly_Laguerre::operator() (double& value) 
{
	double Ln = 1.;
	double product = 1.;
	for (int k = 1; k < base_order; k++) 
	
	{
		product *= base_order - k + 1;
		double ck = factorial(base_order) / (factorial(k) * factorial(base_order - k));
		Ln += std::pow(-1., k) * ck * product * std::pow(value, k);
	
		
	};

	Ln /= factorial(base_order);

	return Ln;
};

int Poly_Laguerre::factorial(int n) 
{
	if (n == 0)
		return 1;
	return n * factorial(n - 1);
};

matrix Poly_Laguerre::operator() (matrix& InputMat) 
{




};