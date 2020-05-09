#pragma once
#include  <vector>
#include  "payoff.hpp"
#include  "sde.hpp"
#include "basisfunction.hpp"


class MonteCarlo
{
	public:
	
		MonteCarlo(){};
		~MonteCarlo(){};
		
		virtual void Simulate(double start, double end, size_t steps)=0;
		
		double GetPrice(double& r, double& T);
		double GetVariance();
		void OptimalNbSimul(const double& errortolerated);
		size_t GetNbSimul();
	
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
		EuropeanBasket_controlvariable(size_t nbSimu, PayOffBasket* Payoff, PayOffBasket* Payoff_control, RandomProcess* diffusion,double inputclosedPrice);
		void Simulate(double start, double end, size_t steps);

	protected:
		PayOffBasket* CPayoff;
		double ExpPriceClsForm;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
class EuropeanBasket_Antithetic: public MonteCarloEuropean
{
	public:

		~EuropeanBasket_Antithetic() {};
		EuropeanBasket_Antithetic(size_t nbSimu, PayOffBasket* Payoff, RandomProcess* diffusion);
		void Simulate(double start, double end, size_t steps);

	protected:
		RandomProcess* x_diffusion;
		matrix paths;
		matrix simulated_price_Anti;
		matrix average_price;
	
};
class EuropeanBasket_Antithetic_CV: public EuropeanBasket_Antithetic
{
	public:

		~EuropeanBasket_Antithetic_CV() {};
		EuropeanBasket_Antithetic_CV(size_t nbSimu, PayOffBasket* Payoff,PayOffBasket* Payoff_control, BSEulerNDAntithetic* diffusion,double inputclosedPrice);
		void Simulate(double start, double end, size_t steps);

	protected:
	
		PayOffBasket* CPayoff;
		double ExpPriceClsForm;
	
};

class AmericanMonteCarlo : public MonteCarlo
{
public:
	AmericanMonteCarlo() {};
	~AmericanMonteCarlo() {};
	AmericanMonteCarlo(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion, 
		std::vector<basis_functions*> polynomial, double df);
	void Simulate(double start, double end, size_t steps);
	matrix GetEarlyExec();

protected:
	std::vector<basis_functions*> Phi;
	double df;
	std::vector<int> early_exec;

};

class AmericanMonteCarlo_controlevariable : public AmericanMonteCarlo
{

public:

	AmericanMonteCarlo_controlevariable() {};
	~AmericanMonteCarlo_controlevariable() {};
	AmericanMonteCarlo_controlevariable(size_t nbSimu, PayOffBasket* InputPayoff, PayOffBasket* CPayoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial, double df);
	void Simulate(double start, double end, size_t steps);

protected:
	PayOffBasket* CPayoff;
	matrix simulated_price_CP;


};

class AmericanMonteCarlo_Antithetic : public AmericanMonteCarlo 
{

public:
	AmericanMonteCarlo_Antithetic() {};
	~AmericanMonteCarlo_Antithetic() {};
	AmericanMonteCarlo_Antithetic(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
		std::vector<basis_functions*> polynomial, double df);
	void Simulate(double start, double end, size_t steps);


protected:
	RandomProcess* x_diffusion;




};