#include "MonteCarlo.h"
#include "Matrix.hpp"

EuropeanVanilla_MonteCarlo::EuropeanVanilla_MonteCarlo(size_t nbSimu, PayOff* Payoff, RandomProcess* diffusion)
{

	MC_price = 0;
	MC_variance = 0.;
	simulated_price.resize(m_Simulation);

};


void EuropeanVanilla_MonteCarlo::Simulate(double start, double end, size_t steps)
{
	for (size_t s = 0; s < m_Simulation; s++)
	{
		m_diffusion->Simulate(start, end, steps);
		matrix paths = m_diffusion->GetAllPaths();
		double maturity_spot = paths(paths.nb_rows() - 1, paths.nb_cols() - 1); //get the spot at maturity

		MC_price = MC_price + Payoff->operator()(maturity_spot) / m_Simulation;
		simulated_price[s] = MC_price;


	}

	for (size_t i = 0; i < simulated_price.size(); i++) 
	{
	
		MC_variance = MC_variance + (simulated_price[i] - MC_price) * (simulated_price[i] - MC_price);
		MC_variance = MC_variance / (simulated_price.size() - 1);
	
	}


};

EuropeanBasket_MonteCarlo::EuropeanBasket_MonteCarlo(size_t nbSimu, PayOffBasket* Payoff, RandomProcess* diffusion)
{

	MC_price = 0;
	MC_variance = 0.;
	simulated_price.resize(m_Simulation);

};

void EuropeanBasket_MonteCarlo::Simulate(double start, double end, size_t steps)
{
	for (size_t s = 0; s < m_Simulation; s++) 
	{
		m_diffusion->Simulate(start, end, steps);
		matrix paths = m_diffusion->GetAllPaths();
		matrix maturity_spot(paths.nb_rows(), 1);

		for (size_t i = 0; i < paths.nb_rows(); i++) 
		{

			maturity_spot(i, 0) = paths(i, paths.nb_cols() - 1);
			//get the vector at maturity 
		}

		MC_price = MC_price + Payoff->operator()(maturity_spot) /m_Simulation;
		simulated_price[s] = MC_price;
	}

	for (size_t i = 0; i < simulated_price.size(); i++) {

		MC_variance = MC_variance + (simulated_price[i] - MC_price) * (simulated_price[i] - MC_price);
		MC_variance = MC_variance / (simulated_price.size() - 1);

	}


};

double EuropeanBasket_MonteCarlo::GetPrice()
{
	return MC_price;
};

double EuropeanBasket_MonteCarlo::GetVariance()
{
	return MC_variance;
};

EuropeanBasket_MonteCarlo_controlvariable::EuropeanBasket_MonteCarlo_controlvariable(size_t nbSimu, 
	PayOffBasket* Payoff, PayOffBasket* Payoff_control, RandomProcess* diffusion):
	EuropeanBasket_MonteCarlo(nbSimu,Payoff,diffusion), CPayoff(Payoff_control)
{

	MC_price = 0;
	MC_variance = 0.;
	simulated_price.resize(m_Simulation);

};

void EuropeanBasket_MonteCarlo::Simulate(double start, double end, size_t steps)
{
	
};
