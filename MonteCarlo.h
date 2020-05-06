#pragma once
#include  <vector>
#include  "payoff.hpp"
#include  "sde.hpp"


class MonteCarlo
{
	public:
	
		MonteCarlo(){};
		~MonteCarlo(){};
		
		virtual void Simulate(double start, double end, size_t steps)=0;
		
		double GetPrice(double& r, double& T);
		double GetVariance();
	
	protected:
	
		size_t m_Simulation;
		PayOffBasket* Payoff;
		RandomProcess* m_diffusion;
		double MC_price;
		double MC_variance;
		matrix simulated_price;

};

class MonteCarloEuropean : public MonteCarlo
{
	public:
	
		MonteCarloEuropean(size_t nbSimu, PayOffBasket* Payoff, RandomProcess* diffusion);
		~MonteCarloEuropean(){};
};
/*
class EuropeanVanilla : public MonteCarloEuropean
{

	public:

		~EuropeanVanilla() {};
		EuropeanVanilla(size_t nbSimu, PayOff* Payoff, RandomProcess* diffusion);
		void Simulate(double start, double end, size_t steps);

};
*/
class EuropeanBasket:public MonteCarloEuropean
{
	
	public:

		~EuropeanBasket() {};
		EuropeanBasket(size_t nbSimu, PayOffBasket* Payoff, RandomProcess* diffusion);
		void Simulate(double start, double end, size_t steps);
};

class EuropeanBasket_controlvariable : public MonteCarloEuropean
{

	public:

		~EuropeanBasket_controlvariable() {};
		EuropeanBasket_controlvariable(size_t nbSimu, PayOffBasket* Payoff, PayOffBasket* Payoff_control, RandomProcess* diffusion);
		void Simulate(double start, double end, size_t steps);

	protected:
		PayOffBasket* CPayoff;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
class EuropeanBasket_Antithetic: public MonteCarloEuropean
{
	public:

		~EuropeanBasket_Antithetic() {};
		EuropeanBasket_Antithetic(size_t nbSimu, PayOffBasket* Payoff, BSEulerNDAntithetic* diffusion);
		void Simulate(double start, double end, size_t steps);

	protected:
		BSEulerNDAntithetic* x_diffusion;
		matrix paths;
	
};