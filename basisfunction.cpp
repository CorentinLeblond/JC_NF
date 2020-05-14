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

	return Ln;
};


Poly_Hermite::Poly_Hermite(int order) :
	basis_functions(order)
{};


double Poly_Hermite::operator() (double& value)
{

	int intpart = base_order/2;

	double Hn = 0.;
	for (int k = 0; k < intpart +1; k++)

	{
		unsigned __int64 facto = (unsigned __int64) factorial(k) * factorial(base_order - 2 * k);
		Hn += std::pow(-1., k) * std::pow(2.*value, base_order- 2*k) /facto;

	}

	Hn = Hn * factorial(base_order);

	return Hn;
};

matrix Poly_Hermite::operator() (matrix& InputMat)
{

	matrix Res(InputMat.nb_rows(), 1);


	for (size_t r = 0; r < InputMat.nb_rows(); r++)
	{

		Res(r, 0) = operator()(InputMat(r, 0));
	}

	return Res;

};


int factorial(int n) 
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
		}
	
		return Res;

};

Laguerre_test::Laguerre_test(int order) :
	basis_functions(order)
{};

double Laguerre_test::operator() (double& value)
{
	double L = 0.;
	switch (base_order) {
	case 0: return 1.;
	case 1: return  1. - value;  
	case 2:  return 0.5 * (2. - 4. * value + value * value);;
	case 3: return (6. - 18. * value + 9. * value * value - value * value * value) / 6.;
	default: return -10000000000000.;
	};
	std::cout << L << std::endl;
};

matrix Laguerre_test::operator() (matrix& InputMat)
{

	matrix Res(InputMat.nb_rows(), 1);


	for (size_t r = 0; r < InputMat.nb_rows(); r++)
	{

		Res(r, 0) = operator()(InputMat(r, 0)); 
	}

	return Res;

};

polynome_simple::polynome_simple(int order) :
	basis_functions(order)
{};

double polynome_simple::operator()(double& value) 
{   
	return std::pow(value, base_order);
};

matrix polynome_simple::operator() (matrix& InputMat)
{

	matrix Res(InputMat.nb_rows(), 1);


	for (size_t r = 0; r < InputMat.nb_rows(); r++)
	{

		Res(r, 0) = operator()(InputMat(r, 0)); 
	}

	return Res;

};