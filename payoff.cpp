
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

double PayOffBasketCall::operator() (matrix inputWeights, matrix inputS,const double& df)  
{
	matrix indexvalue = Weights*S;
	return std::max(indexvalue(0,0) - K*df, 0.0); // Basket Call payoff
};

PayOffControlVarBasketCall::PayOffControlVarBasketCall(matrix inputWeights, matrix inputS,const double& inputK)
{
	Weights = inputWeights;
	S = inputS;
	K = inputK;
};

double PayOffControlVarBasketCall::operator() (matrix Weights, matrix S,const double& df)  
{
	matrix indexvalue = S;
	for(size_t i = 0;i<S.nb_rows();++i)
	{
		indexvalue(i,0) = log(S(i,0));
	}
	indexvalue = Weights*indexvalue;
	return std::max(exp(indexvalue(0,0)) - K*df, 0.0); // Control Variable for Basket Call payoff
};