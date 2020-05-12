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
	double Ln = 0.;
	for (int k = 0; k < base_order+1; k++) 
	
	{
		double ck = factorial(base_order) / (factorial(k) * factorial(base_order - k));
		Ln += std::pow(-1.,k) * ck * std::pow(value, k)/ factorial(k);
	
		
	};

	//Ln /= factorial(base_order);

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

	matrix Res(InputMat.nb_rows(), 1);

	
		for (size_t r = 0; r < InputMat.nb_rows(); r++)
		{
			Res(r, 0) = operator()(InputMat(r, 0)); //
			//std::cout << Res(r, 0) << std::endl;
		}
	
		return Res;

};

Laguerre_test::Laguerre_test(int order) :
	basis_functions(order)
{};

double Laguerre_test::operator() (double& value)
{
	//std::cout << base_order << std::endl;
	double L = 0.;
	switch (base_order) {
	case 0: return 1.;
	case 1: return  1. - value;  
	case 2:  return 0.5 * (2. - 4. * value + value * value);;
	case 3: return (6. - 18. * value + 9. * value * value - value * value * value) / 6.;
	default: return -10000000000000.;
	};
	std::cout << L << std::endl;
	//return L;
};

matrix Laguerre_test::operator() (matrix& InputMat)
{

	matrix Res(InputMat.nb_rows(), 1);


	for (size_t r = 0; r < InputMat.nb_rows(); r++)
	{

		Res(r, 0) = operator()(InputMat(r, 0)); //
		//std::cout << Res(r, 0) << std::endl;
	}

	return Res;

};
