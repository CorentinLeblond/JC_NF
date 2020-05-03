#include "MonteCarlo.h"
#include <iostream>


EuropeanVanilla_MonteCarlo::EuropeanVanilla_MonteCarlo(size_t nbSimu, PayOff* Payoff, 
	RandomProcess* diffusion):
	m_Simulation(nbSimu), Payoff(Payoff), m_diffusion(diffusion)
{

	MC_price = 0;
	MC_variance = 0.;
	simulated_price.Resize(m_Simulation,1);

};


double EuropeanVanilla_MonteCarlo::GetPrice(double& r, double& T)
{
	return exp(-r * T) * MC_price;
};

double EuropeanVanilla_MonteCarlo::GetVariance()
{
	return MC_variance;
};


void EuropeanVanilla_MonteCarlo::Simulate(double start, double end, size_t steps)
{
	for (size_t s = 0; s < m_Simulation; s++)
	{
		m_diffusion->Simulate(start, end, steps);
		matrix paths = m_diffusion->GetAllPaths();
		double maturity_spot = paths(paths.nb_rows() - 1, paths.nb_cols() - 1); //get the spot at maturity
		simulated_price(s,0) = Payoff->operator()(maturity_spot);


	}
	
	MC_price = simulated_price.mean();
	MC_variance = simulated_price.variance();

};

EuropeanBasket_MonteCarlo::EuropeanBasket_MonteCarlo(size_t nbSimu, PayOffBasket* Payoff, 
	RandomProcess* diffusion):
	m_Simulation(nbSimu), Payoff(Payoff), m_diffusion(diffusion)
{

	MC_price = 0;
	MC_variance = 0.;
	simulated_price.Resize(m_Simulation,1);

};

void EuropeanBasket_MonteCarlo::Simulate(double start, double end, size_t steps)
{
	for (size_t s = 0; s < m_Simulation; s++) 
	{
		std::cout << "simulation " << s << std::endl;
		m_diffusion->Simulate(start, end, steps);
		matrix paths = m_diffusion->GetAllPaths();
		matrix maturity_spot(paths.nb_rows(), 1);

		for (size_t i = 0; i < paths.nb_rows(); i++) 
		{

			maturity_spot(i, 0) = paths(i, paths.nb_cols() - 1);
			//get the vector at maturity 
		}

		
		simulated_price(s,0) = Payoff->operator()(maturity_spot);
	}

	MC_price = simulated_price.mean();
	MC_variance = simulated_price.variance();

};

double EuropeanBasket_MonteCarlo::GetPrice(double& r, double& T)
{
	//std::cout << MC_price << std::endl;
	return exp(-r*T)*MC_price;
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
	simulated_price.Resize(m_Simulation,1);

};

void EuropeanBasket_MonteCarlo_controlvariable::Simulate(double start, double end, size_t steps)
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


		simulated_price(s, 0) = Payoff->operator()(maturity_spot)- CPayoff->operator()(maturity_spot);
	}

	MC_price = simulated_price.mean();
	MC_variance = simulated_price.variance();
	
};
