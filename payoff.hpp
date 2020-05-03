
#include <algorithm> 
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <iostream>
#include "Matrix.hpp"

class PayOff 
{
	public:
		//Constructor
		PayOff(){};
		// Virtual destructor to avoid memory leaks when destroying the base and inherited classes
		virtual ~PayOff(){}; 
		// We turn the payoff into a functor, goal is to compute initial condition
		virtual double operator() (const double& S,const double& df = 1) const = 0;
		//virtual double init_cond(const double& S)) const;
		
	protected:
	
		double m_K;
};

//First payoff class created for european call options
//It inherits from Payoff
class PayOffCall : public PayOff 
{
	public:
	
		PayOffCall(const double& K);
		virtual ~PayOffCall() {};

		// Virtual function is now over-ridden (not pure-virtual anymore)
		virtual double operator() (const double& S,const double& df = 1) const;
	  
};
///////////////////////////////////////////////////////////////////////////////////////////

class PayOffBasket
{
	public:
	
		PayOffBasket(){};
		~PayOffBasket(){};
		
		virtual double operator() (matrix S,const double& df = 1) = 0;
		
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

		// Virtual function is now over-ridden (not pure-virtual anymore)
		double operator() (matrix S,const double& df = 1);
	  
};

class PayOffControlVarBasketCall : public PayOffBasket 
{
	public:
	
		PayOffControlVarBasketCall(matrix inputWeights, matrix inputS,const double& inputK);
		~PayOffControlVarBasketCall(){};

		// Virtual function is now over-ridden (not pure-virtual anymore)
		virtual double operator() (matrix S,const double& df = 1);
	  
};