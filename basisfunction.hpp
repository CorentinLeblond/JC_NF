#pragma once 
#include "Matrix.hpp"
class basis_functions 

{

public:

	basis_functions() {};
	~basis_functions() {};
	basis_functions(int order);
	virtual double operator()(double& value) = 0;
	virtual matrix operator()(matrix& InputMat) = 0;

protected:

	int base_order;

};

class Poly_Laguerre : public basis_functions 
{
public:

	Poly_Laguerre() {};
	~Poly_Laguerre() {};
	Poly_Laguerre(int order);
	int factorial(int n);
	double operator()(double& value);
	matrix operator()(matrix& InputMat);

//protected:

};

class Laguerre_test : public basis_functions
{
public:

	Laguerre_test() {};
	~Laguerre_test() {};
	Laguerre_test(int order);
	double operator()(double& value);
	matrix operator()(matrix& InputMat);

};