#include <iostream>
#include "payoff.hpp"
/*
PayOffCall::PayOffCall(const double& K)
{
	m_K = K;
};
	
	//Method to compute the initial condition of the Call option
double PayOffCall::operator() (const double& S,const double& df) const 
{
	return std::max(S-m_K*df, 0.0); // Call payoff
};
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////

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
	//On a la matrice de covariance en paramètre du constructeur, ainsi que ls weights -> on peut calculer la variance du basket
	double variance_basket = 0.;
	double indexvalue = 1.;
	matrix Weights_T(S_.nb_rows(), 1);
	//Spot price + transpose matrix of weights into column to enable compute basket variance later
	for (size_t i = 0; i < S.nb_rows(); ++i)
	{
		//std::cout << "value of S " << indexvalue(i,0) << std::endl;
		indexvalue *= std::pow(S_(i, 0), Weights(0, i));
		Weights_T(i, 0) = Weights(0, i);
		//std::cout << "value of log S " <<  log(S(i, 0)) << std::endl;
	}
	// std::cout << "index val : " << indexvalue << std::endl;
	// std::cout << "Transposed matrix : " << std::endl;
	Weights_T.Print();
	//Variance basket
	Weights *= VarCovar;
	Weights *= Weights_T; //Normalement un scalaire
	// std::cout<<"variance : "<<std::endl;
	// Weights.Print();
	variance_basket = Weights(0, 0);
	// std::cout << "variance : " << variance_basket << std::endl;
	//On peut commencer le pricing
	//variance basket is sigma squared in the bs formula
	double d1 = log(indexvalue / K);
	// std::cout << "d1 : " << d1 << std::endl;
	d1 += (rate + variance_basket / 2) * time;
	// std::cout << "d1 : " << d1 << std::endl;
	d1 /= sqrt(variance_basket * time);
	// std::cout << "d1 : " << d1 << std::endl;
	double d2 = d1 - sqrt(variance_basket * time);
	// std::cout << "d2 : " << d2 << std::endl;
	double N_d1 = normalCDF(d1);
	double N_d2 = normalCDF(d2);
	// std::cout << "N_d1 : " << N_d1 << std::endl;
	// std::cout << "N_d2 : " << N_d2 << std::endl;
	return (indexvalue * N_d1 - df * K * N_d2);
};

double normalCDF(double x)
{
	return std::erfc(-x / std::sqrt(2)) / 2;
};