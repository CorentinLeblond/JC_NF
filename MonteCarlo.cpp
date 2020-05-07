#include "MonteCarlo.h"
#include <iostream>

//MC
//////////////////////////////////////////////////////////

double MonteCarlo::GetPrice(double& r, double& T)
{
	return exp(-r * T) * MC_price;
};

double MonteCarlo::GetVariance()
{
	return MC_variance;
};
void MonteCarlo::OptimalNbSimul(const double& errortolerated)
{
	//The two-sided 99% quantile is:
	double quantile = 2.57;
	
	//It returns the minimum number of simulations required to make the estimated price below the error tolerated with confidence level 99%
	m_Simulation = (size_t) ((quantile*quantile*MC_variance)/(errortolerated*errortolerated)); // Ptetre changer le type de la méthode comme on veut un nombre de simulation
	
	//The process is tho run a first companion MC on a few simulations, then estimate the variance, then definie minimul number and run the MC over this computed number
};
//MonteCarlo European
/////////////////////////////////////////////////////////
MonteCarloEuropean::MonteCarloEuropean(size_t nbSimu, PayOffBasket* inputPayoff, RandomProcess* diffusion)
{
	Payoff = inputPayoff;
	m_Simulation = nbSimu;
	m_diffusion = diffusion;
	MC_price = 0;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu,1);
};

//Single asset Vanilla European Options
////////////////////////////////////////////////////////

/*
EuropeanVanilla::EuropeanVanilla(size_t nbSimu, PayOff* inputPayoff, RandomProcess* diffusion):
	MonteCarloEuropean(nbSimu,inputPayoff,diffusion)
{
};

void EuropeanVanilla::Simulate(double start, double end, size_t steps)
{
	// std::cout << "MC European" << std::endl;
	for (size_t s = 0; s < m_Simulation; s++)
	{
		//std::cout << "simulation " << s << std::endl;
		m_diffusion->Simulate(start, end, steps);
		matrix paths = m_diffusion->GetAllPaths();
		double maturity_spot = paths(paths.nb_rows() - 1, paths.nb_cols() - 1); //get the spot at maturity
		simulated_price(s,0) = Payoff->operator()(maturity_spot);

	}
	
	MC_price = simulated_price.mean();
	MC_variance = simulated_price.variance();

};
*/
//European Basket Options
/////////////////////////////////////////////////////////////

EuropeanBasket::EuropeanBasket(size_t nbSimu, PayOffBasket* inputPayoff, RandomProcess* diffusion):
	MonteCarloEuropean(nbSimu,inputPayoff,diffusion)
{
};

void EuropeanBasket::Simulate(double start, double end, size_t steps)
{

	for(size_t s = 0; s < m_Simulation; ++s) 
	{

		m_diffusion->Simulate(start, end, steps);
		
		matrix paths = m_diffusion->GetAllPaths();
		
		matrix maturity_spot(paths.nb_rows(), 1);
		
		std::cout << "MC mat spot"<< std::endl;
		
		for (size_t i = 0; i < paths.nb_rows(); i++) 
		{
			//std::cout << paths(i, paths.nb_cols() - 1) << std::endl;
			maturity_spot(i, 0) = paths(i, paths.nb_cols() - 1);

			//get the vector at maturity 
		}

		simulated_price(s,0) = Payoff->operator()(maturity_spot);
	}

	MC_price = simulated_price.mean();
	MC_variance = simulated_price.variance();

};

EuropeanBasket_controlvariable::EuropeanBasket_controlvariable(size_t nbSimu, PayOffBasket* inputPayoff, PayOffBasket* Payoff_control, RandomProcess* diffusion):
	MonteCarloEuropean(nbSimu,inputPayoff,diffusion),
	CPayoff(Payoff_control)
{
};

void EuropeanBasket_controlvariable::Simulate(double start, double end, size_t steps)
{
	std::cout << "MC European Basket and CV" << std::endl;
	matrix CVprice(m_Simulation, 1);

	for (size_t s = 0; s < m_Simulation; s++)
	{
		//std::cout << "simulation " << s << std::endl;
		m_diffusion->Simulate(start, end, steps);
		matrix paths = m_diffusion->GetAllPaths();
		//std::cout << "Path matrix is  " << s << std::endl;
		//paths.Print();
		matrix maturity_spot(paths.nb_rows(), 1);

		for (size_t i = 0; i < paths.nb_rows(); i++)
		{

			maturity_spot(i, 0) = paths(i, paths.nb_cols() - 1);
			//get the vector at maturity 
		}

		//std::cout << "mat spot at " << s << std::endl;
		//maturity_spot.Print();

		simulated_price(s, 0) = Payoff->operator()(maturity_spot)- CPayoff->operator()(maturity_spot);
		CVprice(s, 0) = CPayoff->operator()(maturity_spot);


		//if (Payoff->operator()(maturity_spot) > CPayoff->operator()(maturity_spot)) 
		//{std::cout << "simu " << s << "difference positive" << std::endl;}
		

	}

	
	MC_price = simulated_price.mean() + CVprice.mean();
	MC_variance = simulated_price.variance();
	
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////
EuropeanBasket_Antithetic::EuropeanBasket_Antithetic(size_t nbSimu, PayOffBasket* inputPayoff, BSEulerNDAntithetic* diffusion):
	MonteCarloEuropean(nbSimu,inputPayoff,diffusion),
	x_diffusion(diffusion)
{
};

void EuropeanBasket_Antithetic::Simulate(double start, double end, size_t steps)
{
	std::cout << "MC European Basket"<< std::endl;

	for (size_t s = 0; s < m_Simulation; s++) 
	{
		//std::cout << "simulation " << s << std::endl;
		if(s % 2 ==0)
		{
			x_diffusion->Simulate(start, end, steps);
			paths = x_diffusion->GetAllPaths();
		}
		else
		{
			paths = x_diffusion->GetAllPathsAnti();
		}
				
		//std::cout << "path" << std::endl;
		//paths.Print();
		//std::cout << paths.nb_rows() << std::endl;
		matrix maturity_spot(paths.nb_rows(), 1);

		for (size_t i = 0; i < paths.nb_rows(); i++) 
		{
			//std::cout << paths(i, paths.nb_cols() - 1) << std::endl;
			maturity_spot(i, 0) = paths(i, paths.nb_cols() - 1);
			//std::cout << "index" << std::endl;
			
			//get the vector at maturity 
		}

		//maturity_spot.Print();
		
		simulated_price(s,0) = Payoff->operator()(maturity_spot);
	}

	MC_price = simulated_price.mean();
	MC_variance = simulated_price.variance();

};