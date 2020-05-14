#pragma once
#include <algorithm> 
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Matrix.hpp"


class PayOffBasket
{
	public:
	
		PayOffBasket(){};
		~PayOffBasket(){};
		matrix GetWeights();
		
		virtual double operator() (matrix S,const double& df = 1) = 0;
		virtual double operator() (double I, const double& df = 1) = 0;
		
	protected:
	
		matrix Weights;
		matrix S;
		double K;
};
class PayOffBasketCall : public PayOffBasket 
{
	public:
	
		PayOffBasketCall(matrix inputWeights, matrix inputS,const double& inputK);
		~PayOffBasketCall() {};

		double operator() (matrix S,const double& df = 1);
		double operator() (double I, const double& df = 1);
	  
};

class PayOffControlVarBasketCall : public PayOffBasket 
{
	public:
	
		PayOffControlVarBasketCall(matrix inputWeights, matrix inputS,const double& inputK);
		~PayOffControlVarBasketCall(){};

		virtual double operator() (matrix S,const double& df = 1);
		virtual double operator() (double I, const double& df = 1);
	  
};

class ClosedFormulaBasketCall : public PayOffBasket
{
public:

	ClosedFormulaBasketCall(matrix inputWeights, matrix inputS, matrix inputVarCovar, const double& inputK, const double& inputrate, const double& inputtime);
	~ClosedFormulaBasketCall() {};

	double operator() (matrix S, const double& df = 1);
	double operator() (double I, const double& df = 1);
	
private:

	matrix VarCovar;
	double rate;
	double time;
};
double normalCDF(double x);
