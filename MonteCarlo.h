#pragma once
#include  <vector>
#include  "payoff.hpp"
#include  "sde.hpp"


class EuropeanVanilla_MonteCarlo
{

public:

	EuropeanVanilla_MonteCarlo() {};
	~EuropeanVanilla_MonteCarlo() {};
	EuropeanVanilla_MonteCarlo(size_t nbSimu, PayOff* Payoff, RandomProcess* diffusion);
	void Simulate(double start, double end, size_t steps);
	double GetPrice(double&r, double&T);
	double GetVariance();

protected:
	size_t m_Simulation;
	PayOff* Payoff;
	RandomProcess* m_diffusion;
	double MC_price;
	double MC_variance;
	matrix simulated_price;


};

class EuropeanBasket_MonteCarlo 
{
	
public:

	EuropeanBasket_MonteCarlo() {};
	~EuropeanBasket_MonteCarlo() {};
	EuropeanBasket_MonteCarlo(size_t nbSimu, PayOffBasket* Payoff, RandomProcess* diffusion);
	void Simulate(double start, double end, size_t steps);
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

class EuropeanBasket_MonteCarlo_controlvariable : public EuropeanBasket_MonteCarlo
{

public:

	EuropeanBasket_MonteCarlo_controlvariable() {};
	~EuropeanBasket_MonteCarlo_controlvariable() {};
	EuropeanBasket_MonteCarlo_controlvariable(size_t nbSimu, PayOffBasket* Payoff, PayOffBasket* Payoff_control, RandomProcess* diffusion);
	void Simulate(double start, double end, size_t steps);

protected:
	size_t m_Simulation;
	PayOffBasket* Payoff;
	PayOffBasket* CPayoff;
	RandomProcess* m_diffusion;

};
