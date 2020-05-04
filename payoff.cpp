#include <iostream>
#include "payoff.hpp"

PayOffCall::PayOffCall(const double& K)
{
	m_K = K;
};
	
	//Method to compute the initial condition of the Call option
double PayOffCall::operator() (const double& S,const double& df) const 
{
	return std::max(S-m_K*df, 0.0); // Call payoff
};

/////////////////////////////////////////////////////////////////////////////////////////////////////

PayOffBasketCall::PayOffBasketCall(matrix inputWeights, matrix inputS,const double& inputK)
{
	Weights = inputWeights;
	S = inputS;
	K = inputK;
};

double PayOffBasketCall::operator() (matrix inputS,const double& df)  
{
	matrix indexvalue = Weights*S;

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