#include <iostream>
#include "payoff.hpp"


/* This file contains all payoff related algorithm, for the basket call and control variate*/
matrix PayOffBasket::GetWeights() { return Weights; };

PayOffBasketCall::PayOffBasketCall(matrix inputWeights, matrix inputS,const double& inputK)
{
	Weights = inputWeights;
	S = inputS;
	K = inputK;
};

double PayOffBasketCall::operator() (matrix inputS,const double& df)  
{
	matrix indexvalue = Weights*inputS;

	//inputS.Print();
	//std::cout << "index bskt " << std::endl;
	//indexvalue.Print();

	return std::max(indexvalue(0,0) - K*df, 0.0); // Basket Call payoff
};

double PayOffBasketCall::operator() (double I, const double& df)
{

	return std::max(I - K * df, 0.0);


};


PayOffControlVarBasketCall::PayOffControlVarBasketCall(matrix inputWeights, matrix inputS,const double& inputK)
{
	Weights = inputWeights;
	S = inputS;
	K = inputK;
};

double PayOffControlVarBasketCall::operator() (matrix S,const double& df)  
{
	matrix indexvalue = S;
	for(size_t i = 0;i<S.nb_rows();++i)
	{
		//std::cout << "value of S " << indexvalue(i,0) << std::endl;
		indexvalue(i,0) = log(S(i,0));
		//std::cout << "value of log S " <<  log(S(i, 0)) << std::endl;
	}
	indexvalue = Weights*indexvalue;

	//std::cout << "index CV " << std::endl;
	//indexvalue.Print();
	return std::max(exp(indexvalue(0,0)) - K*df, 0.0); // Control Variable for Basket Call payoff
};

double  PayOffControlVarBasketCall::operator() (double I, const double& df)
{

	return std::max(exp(I) - K * df, 0.0);


};


ClosedFormulaBasketCall::ClosedFormulaBasketCall(matrix inputWeights, matrix inputS, matrix inputVarCovar, const double& inputK,
	const double& inputrate, const double& inputtime) :
	VarCovar(inputVarCovar),
	rate(inputrate),
	time(inputtime)
{
	Weights = inputWeights;
	S = inputS;
	K = inputK;
};
double ClosedFormulaBasketCall::operator() (matrix S_, const double& df)
{
	double variance_basket = 0.;
	double squared_var = 0.;
	double indexvalue = 1.;
	matrix Weights_T(S_.nb_rows(), 1);
	Weights_T = transpose(Weights);
	//Spot price + transpose matrix of weights into column to enable compute basket variance later
	for (size_t i = 0; i < S.nb_rows(); ++i)
	{
		indexvalue *= std::pow(S_(i, 0), Weights(0, i));
		squared_var += VarCovar(i, i)* Weights(0, i);
	}

	Weights *= VarCovar;
	Weights *= Weights_T; //Normalement un scalaire

	variance_basket = Weights(0, 0);

	double d1 = log(indexvalue / K);
	d1 += (rate + 0.5*variance_basket -0.5*squared_var / 2) * time;

	d1 /= sqrt(variance_basket * time);

	double d2 = d1 - sqrt(variance_basket * time);

	double N_d1 = normalCDF(d1);
	double N_d2 = normalCDF(d2);

	return (indexvalue * N_d1 - df * K * N_d2);
};


double ClosedFormulaBasketCall::operator() (double I, const double& df)
{

	return std::max(I - K * df, 0.0);

};

double normalCDF(double x)
{
	return std::erfc(-x / std::sqrt(2)) / 2;
};