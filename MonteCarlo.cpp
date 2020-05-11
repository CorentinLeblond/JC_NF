#include "MonteCarlo.h"
#include <iostream>

//MC
//////////////////////////////////////////////////////////

double MonteCarlo::GetPrice()
{
	return MC_price;
};

size_t MonteCarlo::GetNbSimul()
{
	return m_Simulation;
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
	r = diffusion->Get_rate();
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
	std::cout << "MC European Basket" << std::endl;
	for(size_t s = 0; s < m_Simulation; ++s) 
	{

		m_diffusion->Simulate(start, end, steps);
		
		matrix paths = m_diffusion->GetAllPaths();
		
		matrix maturity_spot(paths.nb_rows(), 1);
		
		//std::cout << "MC mat spot"<< std::endl;
		
		for (size_t i = 0; i < paths.nb_rows(); i++) 
		{
			//std::cout << paths(i, paths.nb_cols() - 1) << std::endl;
			maturity_spot(i, 0) = paths(i, paths.nb_cols() - 1);

			//get the vector at maturity 
		}

		simulated_price(s,0) = Payoff->operator()(maturity_spot);
	}

	MC_price = exp(-r*end)*simulated_price.mean();
	MC_variance = simulated_price.variance();

};

EuropeanBasket_controlvariable::EuropeanBasket_controlvariable(size_t nbSimu, PayOffBasket* inputPayoff, PayOffBasket* Payoff_control,
																RandomProcess* diffusion,double inputclosedPrice)
	:
	MonteCarloEuropean(nbSimu,inputPayoff,diffusion),
	CPayoff(Payoff_control),
	ExpPriceClsForm(inputclosedPrice)
{
};

void EuropeanBasket_controlvariable::Simulate(double start, double end, size_t steps)
{
	std::cout << "MC European Basket and CV" << std::endl;
	// matrix CVprice(m_Simulation, 1);

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

		simulated_price(s, 0) = Payoff->operator()(maturity_spot)- CPayoff->operator()(maturity_spot)+ExpPriceClsForm;
		// CVprice(s, 0) = CPayoff->operator()(maturity_spot);


		//if (Payoff->operator()(maturity_spot) > CPayoff->operator()(maturity_spot)) 
		//{std::cout << "simu " << s << "difference positive" << std::endl;}
	}

	
	// MC_price = simulated_price.mean() + CVprice.mean();
	MC_price = exp(-r*end)*simulated_price.mean();
	MC_variance = simulated_price.variance();
	
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////
EuropeanBasket_Antithetic::EuropeanBasket_Antithetic(size_t nbSimu, PayOffBasket* inputPayoff, RandomProcess* diffusion):
	MonteCarloEuropean(nbSimu,inputPayoff,diffusion),
	x_diffusion(diffusion)
{
	simulated_price_Anti.Resize(m_Simulation/2,1);
	average_price.Resize(m_Simulation/2,1);
};

void EuropeanBasket_Antithetic::Simulate(double start, double end, size_t steps)
{
	std::cout << "MC European Basket and Anti"<< std::endl;
	
	simulated_price.Resize(m_Simulation/2,1);
	// simulated_price_Anti.Resize(m_Simulation/2,1);
	size_t k = 0;
	
	std::cout << "MC European Basket antithetic"<< std::endl;

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
		if(s % 2 ==0)
		{
			simulated_price(k,0) = Payoff->operator()(maturity_spot);
		}
		else
		{
			simulated_price_Anti(k,0) = Payoff->operator()(maturity_spot);
			average_price(k,0) = 0.5*(simulated_price_Anti(k,0) + simulated_price(k,0));
			k += 1;
		}		
		
	}

	MC_price = exp(-r*end)*average_price.mean();
	MC_variance = average_price.variance();

};

EuropeanBasket_Antithetic_CV::EuropeanBasket_Antithetic_CV(size_t nbSimu, PayOffBasket* Payoff,PayOffBasket* Payoff_control,
														BSEulerNDAntithetic* diffusion,double inputclosedPrice)
	:
	EuropeanBasket_Antithetic(nbSimu,Payoff,diffusion),
	CPayoff(Payoff_control),
	ExpPriceClsForm(inputclosedPrice)
{
};
//////////////////////////////////////////////////////////////////////////////////////////////////
void EuropeanBasket_Antithetic_CV::Simulate(double start, double end, size_t steps)
{
	// std::cout << "MC European Basket and CV" << std::endl;
	std::cout << "MC European Basket and Anti with CV"<< std::endl;
	simulated_price.Resize(m_Simulation/2,1);
	// simulated_price_Anti.Resize(m_Simulation/2,1);
	size_t k = 0;
	
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
		if(s % 2 ==0)
		{
			simulated_price(k,0) = Payoff->operator()(maturity_spot)-CPayoff->operator()(maturity_spot)+ExpPriceClsForm;
		}
		else
		{
			simulated_price_Anti(k,0) = Payoff->operator()(maturity_spot)- CPayoff->operator()(maturity_spot)+ExpPriceClsForm;
			average_price(k,0) = 0.5*(simulated_price_Anti(k,0) + simulated_price(k,0));
			k += 1;
		}		
		
	}

	MC_price = exp(-r*end)*average_price.mean();
	MC_variance = average_price.variance();
	
};	

/////////////////// US ///////////////////////////////////////////////////////


AmericanMonteCarlo::AmericanMonteCarlo(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial) :
	Phi(polynomial)
{
	r = diffusion->Get_rate();
	Payoff = InputPayoff;
	m_Simulation = nbSimu;
	m_diffusion = diffusion;
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu,1);
};

matrix AmericanMonteCarlo::C_Hat_regression(matrix Index_time_t, matrix Value)

{
	//Index is a column vector with rows of number of simulations done and value should have the same size !!! 
	size_t size_simu = Index_time_t.nb_rows();
	matrix Phit(size_simu, Phi.size());
	matrix poly_order_n(size_simu, 1);
	matrix Phit_T(Phi.size(), size_simu);
	matrix SDP_Phi(Phi.size(), Phi.size());
	matrix Phi_V(Phi.size(), 1);
	matrix inv_SDP_Phi(Phi.size(), Phi.size());
	matrix Beta_hat(Phi.size(), 1);
	matrix C_hat(size_simu, 1);

	for (size_t b = 0; b < Phi.size(); b++)

	{
		poly_order_n = Phi[b]->operator()(Index_time_t);

		for (size_t i = 0; i < poly_order_n.nb_rows(); i++)
		{
			Phit(i, b) = poly_order_n(i, 0);
		}
	}

	Phit_T = transpose(Phit);

	SDP_Phi = Phit_T * Phit;

	Phi_V = Phit_T * Value;

	inv_SDP_Phi = Inverse(SDP_Phi, SDP_Phi.nb_cols());

	//std::cout << "Inverse of Phi " << std::endl;
	//inv_SDP_Phi.Print();

	Beta_hat = inv_SDP_Phi * Phi_V;

	C_hat = Phit * Beta_hat;

	Phit.Clear();
	poly_order_n.Clear();
	Phit_T.Clear();
	SDP_Phi.Clear();
	Phi_V.Clear();
	inv_SDP_Phi.Clear();
	Beta_hat.Clear();

	//C_hat.Print();
	
	return C_hat;

};

AmericanMonteCarlo_basket::AmericanMonteCarlo_basket(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial) :
	AmericanMonteCarlo(nbSimu,InputPayoff,diffusion,polynomial)
{
	r = diffusion->Get_rate();
	Payoff = InputPayoff;
	m_Simulation = nbSimu;
	m_diffusion = diffusion;
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu, 1);
};



void AmericanMonteCarlo_basket::Simulate(double start, double end, size_t steps)
{

	
	std::cout << "American MC LS" << std::endl;
	simulated_price.Resize(m_Simulation, 1);
	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation, steps);
	matrix ITM(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix It(m_Simulation,1);
	matrix Phit(m_Simulation, Phi.size());
	matrix poly_order_n(m_Simulation, 1);
	matrix Phit_T(Phi.size(),m_Simulation);
	matrix SDP_Phi(Phi.size(),Phi.size());
	matrix Phi_V(Phi.size(),1);
	matrix inv_SDP_Phi(Phi.size(), Phi.size());
	matrix Beta_hat(Phi.size(),1);
	matrix C_hat(m_Simulation, 1);
	matrix chm(1,steps);
	matrix V(m_Simulation, 1);
	matrix paths;
	

	for (size_t s = 0; s < m_Simulation; s++)
	{
		m_diffusion->Simulate(start, end, steps);
		paths = m_diffusion->GetAllPaths();
		
		chm = weights * paths;

		for (size_t i = 0; i < steps; i++) 
		{
			Index(s, i) = chm(0,i);
		
		}
		
		simulated_price(s, 0) = Payoff->operator()(Index(s,Index.nb_cols()-1));
		//populate the vector at maturity of the payoff

	}

	double df = exp(-r * m_diffusion->Get_Dt());

	// Regression part 
	for(size_t t = steps - 2; t>0; t--) 
	{
	
		
		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			V(i, 0) = df * simulated_price(i, 0);
			It(i, 0) = Index(i, t); //take index at t-1;

			if (Payoff->operator()(It(i,0)) > 0)
			{ ITM(i, 0) = 1;  }
			else { ITM(i, 0) = -1; };
		}

		for (size_t b = 0; b < Phi.size(); b++) 
		
		{	
			poly_order_n = Phi[b]->operator()(It);

			for (size_t i = 0; i < poly_order_n.nb_rows(); i++)
			{
				Phit(i, b) = poly_order_n(i,0);
			}
			

		}

		Phit_T = transpose(Phit);

		SDP_Phi = Phit_T * Phit;

		Phi_V = Phit_T * V;

		inv_SDP_Phi = Inverse(SDP_Phi,SDP_Phi.nb_cols());

		Beta_hat = inv_SDP_Phi * Phi_V;

		C_hat = Phit * Beta_hat;

		for (size_t i = 0; i < C_hat.nb_rows(); i++) 
		{
		
			if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(i, t)) > C_hat(i, 0)))
			{
				simulated_price(i,0) = Payoff->operator()(Index(i, t));
				early_exec.push_back(i);
			}
			
			{
			
			
				simulated_price(i, 0) = df * simulated_price(i, 0);
			}
		
		
		}
	}

	MC_price = df* simulated_price.mean();

	MC_variance = simulated_price.variance();

};

matrix AmericanMonteCarlo::GetEarlyExec() 
{
	matrix res(early_exec.size(), 1);

	for (size_t i = 0; i < early_exec.size(); i++)
	{

		res(i, 0) = early_exec[i];

	};

	return res;

};


AmericanMonteCarlo_basket_controlevariable::AmericanMonteCarlo_basket_controlevariable(size_t nbSimu, PayOffBasket* InputPayoff, PayOffBasket* CPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial, double closed_form_price)
	:AmericanMonteCarlo(nbSimu,InputPayoff,diffusion,polynomial), CPayoff(CPayoff), ExpPriceClsForm(closed_form_price)
{
	r = diffusion->Get_rate();
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu, 1);
	simulated_price_CP.Resize(nbSimu,1); //to compute the mean of the control variate payoff 
};


void AmericanMonteCarlo_basket_controlevariable::Simulate(double start, double end, size_t steps)
{

	std::cout << "American MC LS with Control Variable" << std::endl;

	simulated_price.Resize(m_Simulation, 1);
	simulated_price_CP.Resize(m_Simulation, 1);
	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation, steps);
	matrix ITM(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix It(m_Simulation, 1);
	matrix Phit(m_Simulation, Phi.size());
	matrix poly_order_n(m_Simulation, 1);
	matrix Phit_T(Phi.size(), m_Simulation);
	matrix SDP_Phi(Phi.size(), Phi.size());
	matrix Phi_V(Phi.size(), 1);
	matrix inv_SDP_Phi(Phi.size(), Phi.size());
	matrix Beta_hat(Phi.size(), 1);
	matrix C_hat(m_Simulation, 1);
	matrix chm(1, steps);
	matrix V(m_Simulation, 1);
	matrix paths;


	matrix log_Index(m_Simulation, steps);
	matrix ITM2(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix It2(m_Simulation, 1);
	matrix Phit_2(m_Simulation, Phi.size());
	matrix poly_order_n_2(m_Simulation, 1);
	matrix Phit_T_2(Phi.size(), m_Simulation);
	matrix SDP_Phi_2(Phi.size(), Phi.size());
	matrix Phi_V_2(Phi.size(), 1);
	matrix inv_SDP_Phi_2(Phi.size(), Phi.size());
	matrix Beta_hat_2(Phi.size(), 1);
	matrix C_hat_2(m_Simulation, 1);
	matrix chm_2(1, steps);
	matrix V_2(m_Simulation, 1);
	matrix log_spot;

	matrix CV(m_Simulation, 1);

	for (size_t s = 0; s < m_Simulation; s++)
	{
		m_diffusion->Simulate(start, end, steps);
		paths = m_diffusion->GetAllPaths();
		log_spot = paths;

		for (size_t i = 0; i < paths.nb_rows(); i++)
		{
			for (size_t j = 0; j < paths.nb_cols(); j++)
			{
				log_spot(i, j) = log(paths(i, j));
			}
		}
		chm = weights * paths;
		chm_2 = weights * log_spot;

		for (size_t i = 0; i < steps; i++)
		{
			Index(s, i) = chm(0, i);
			log_Index(s, i) = chm_2(0, i);
		}

		simulated_price(s, 0) = Payoff->operator()(Index(s, Index.nb_cols() - 1));
		simulated_price_CP(s,0) = CPayoff->operator()(log_Index(s, Index.nb_cols() - 1));
		//populate the vector at maturity of the payoff

	}

	double df = exp(-r * m_diffusion->Get_Dt());
	// Regression part 

	for (size_t t = steps - 2; t > 0; t--)
	{
		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			V(i, 0) = df * simulated_price(i, 0);
			V_2(i, 0) = df * simulated_price_CP(i, 0);
			It(i, 0) = Index(i, t); //take index at t-1;
			It2(i, 0) = log_Index(i, t);
			if (Payoff->operator()(It(i, 0)) > 0)
			{
				ITM(i, 0) = 1;
			}
			else { ITM(i, 0) = -1; };

			if (CPayoff->operator()(It2(i, 0)) > 0)
			{
				ITM2(i, 0) = 1;
			}
			else { ITM2(i, 0) = -1; };



		}

		for (size_t b = 0; b < Phi.size(); b++)

		{
			poly_order_n = Phi[b]->operator()(It);
			poly_order_n_2 = Phi[b]->operator()(It2);

			for (size_t i = 0; i < poly_order_n.nb_rows(); i++)
			{
				Phit(i, b) = poly_order_n(i, 0);
				Phit_2(i, b) = poly_order_n_2(i, 0);
			}


		}

		Phit_T = transpose(Phit);
		Phit_T_2 = transpose(Phit_2);

		SDP_Phi = Phit_T * Phit;
		SDP_Phi_2 = Phit_T_2 * Phit_2;

		Phi_V = Phit_T * V;
		Phi_V_2 = Phit_T_2 * V_2;

		inv_SDP_Phi = Inverse(SDP_Phi, SDP_Phi.nb_cols());
		inv_SDP_Phi_2 = Inverse(SDP_Phi_2, SDP_Phi_2.nb_cols());

		Beta_hat = inv_SDP_Phi * Phi_V;
		Beta_hat_2 = inv_SDP_Phi_2 * Phi_V_2;


		C_hat = Phit * Beta_hat;
		C_hat_2 = Phit_2 * Beta_hat_2;

		for (size_t i = 0; i < C_hat.nb_rows(); i++)
		{

			if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(i, t)) > C_hat(i, 0)))
			{
				simulated_price(i, 0) = Payoff->operator()(Index(i, t));
			}

			{


				simulated_price(i, 0) = df * simulated_price(i, 0);
			}




			if ((ITM2(i, 0) == 1) && (CPayoff->operator()(log_Index(i, t)) > C_hat_2(i, 0)))
			{
				simulated_price_CP(i, 0) = CPayoff->operator()(log_Index(i, t));
			}

			{


				simulated_price_CP(i, 0) = df * simulated_price_CP(i, 0);
			}

			CV(i, 0) = simulated_price(i, 0) - simulated_price_CP(i, 0) + ExpPriceClsForm;;


		}
	} 


	MC_price = df* CV.mean();

	MC_variance = CV.variance();

};

AmericanMonteCarlo_basket_Antithetic::AmericanMonteCarlo_basket_Antithetic(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial) 
	:AmericanMonteCarlo(nbSimu,InputPayoff,diffusion,polynomial), x_diffusion(diffusion)
{
	r = diffusion->Get_rate();
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu/2, 1);
	simulated_price_anti.Resize(nbSimu/2, 1);
	average_price.Resize(nbSimu/2, 1);

};

void AmericanMonteCarlo_basket_Antithetic::Simulate(double start, double end, size_t steps) 
{

	std::cout << "American MC with antithetic" << std::endl;


	simulated_price.Resize(m_Simulation/2, 1);
	simulated_price_anti.Resize(m_Simulation/2, 1);
	average_price.Resize(m_Simulation/2, 1);

	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation/2, steps);
	matrix ITM(m_Simulation/2, 1); //need to separate which paths are ITM at each iteration 
	matrix It(m_Simulation/2, 1);
	matrix Phit(m_Simulation/2, Phi.size());
	matrix poly_order_n(m_Simulation/2, 1);
	matrix Phit_T(Phi.size(), m_Simulation/2);
	matrix SDP_Phi(Phi.size(), Phi.size());
	matrix Phi_V(Phi.size(), 1);
	matrix inv_SDP_Phi(Phi.size(), Phi.size());
	matrix Beta_hat(Phi.size(), 1);
	matrix C_hat(m_Simulation/2, 1);
	matrix chm(1, steps);
	matrix V(m_Simulation/2, 1);
	matrix paths;
	matrix paths_anti;


	matrix chm_anti(1, steps);
	matrix Index_anti(m_Simulation/2, steps);
	matrix ITM_anti(m_Simulation/2, 1); //need to separate which paths are ITM at each iteration 
	matrix It_anti(m_Simulation/2, 1);
	matrix Phit_anti(m_Simulation/2, Phi.size());
	matrix poly_order_n_anti(m_Simulation/2, 1);
	matrix Phit_T_anti(Phi.size(), m_Simulation/2);
	matrix SDP_Phi_anti(Phi.size(), Phi.size());
	matrix Phi_V_anti(Phi.size(), 1);
	matrix inv_SDP_Phi_anti(Phi.size(), Phi.size());
	matrix Beta_hat_anti(Phi.size(), 1);
	matrix C_hat_anti(m_Simulation/2, 1);
	matrix V_anti(m_Simulation / 2, 1);

	size_t k = 0;
	size_t k_a = 0;

	for (size_t s = 0; s < m_Simulation; s++)
	{

		//std::cout << "simulation " << s << std::endl;
		if (s % 2 == 0)
		{
			x_diffusion->Simulate(start, end, steps);
			paths = x_diffusion->GetAllPaths();
			chm = weights * paths;

			for (size_t i = 0; i < steps; i++)
			{
				Index(k, i) = chm(0, i);

			}

			simulated_price(k, 0) = Payoff->operator()(Index(k, Index.nb_cols() - 1));

			k += 1;

			

		}
		else
		{
			paths_anti = x_diffusion->GetAllPathsAnti();
			chm_anti = weights * paths_anti;


			for (size_t i = 0; i < steps; i++)
			{
				Index_anti(k_a, i) = chm_anti(0, i);

			}

			simulated_price_anti(k_a, 0) = Payoff->operator()(Index_anti(k_a, Index_anti.nb_cols() - 1));
			k_a += 1;
		}
	}

	double df = exp(-r * m_diffusion->Get_Dt());
	// Regression part 

	for (size_t t = steps - 2; t > 0; t--)
	{

		//std::cout << "regression " << t << std::endl;
		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			V(i, 0) = df * simulated_price(i, 0);
			It(i, 0) = Index(i, t); //take index at t-1;

			if (Payoff->operator()(It(i, 0)) > 0)
			{
				ITM(i, 0) = 1;
			}
			else { ITM(i, 0) = -1; };

			V_anti(i, 0) = df * simulated_price_anti(i, 0);
			It_anti(i, 0) = Index_anti(i, t); //take index at t-1;

			if (Payoff->operator()(It_anti(i, 0)) > 0)
			{
				ITM_anti(i, 0) = 1;
			}
			else { ITM_anti(i, 0) = -1; };
		}

		for (size_t b = 0; b < Phi.size(); b++)

		{
			poly_order_n = Phi[b]->operator()(It);
			poly_order_n_anti = Phi[b]->operator()(It_anti);

			for (size_t i = 0; i < poly_order_n.nb_rows(); i++)
			{
				Phit(i, b) = poly_order_n(i, 0);

				Phit_anti(i, b) = poly_order_n_anti(i, 0);
			}


		}

		Phit_T = transpose(Phit);
		Phit_T_anti = transpose(Phit_anti);

		SDP_Phi = Phit_T * Phit;

		SDP_Phi_anti = Phit_T_anti * Phit_anti;

		Phi_V = Phit_T * V;
		Phi_V_anti = Phit_T_anti * V_anti;

		//std::cout << "simulation Phi" << std::endl;
		inv_SDP_Phi = Inverse(SDP_Phi, SDP_Phi.nb_cols());
		//inv_SDP_Phi.Print();
		//std::cout << "simulation Phi anti" << std::endl;
		inv_SDP_Phi_anti = Inverse(SDP_Phi_anti, SDP_Phi_anti.nb_cols());
		//inv_SDP_Phi_anti.Print();

		Beta_hat = inv_SDP_Phi * Phi_V;
		Beta_hat_anti = inv_SDP_Phi_anti * Phi_V_anti;

		C_hat = Phit * Beta_hat;
		C_hat_anti = Phit_anti * Beta_hat_anti;

		for (size_t i = 0; i < C_hat.nb_rows(); i++)
		{

			if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(i, t)) > C_hat(i, 0)))
			{
				simulated_price(i, 0) = Payoff->operator()(Index(i, t));
				//early_exec.push_back(i);
			}

			{


				simulated_price(i, 0) = df * simulated_price(i, 0);
			}


			if ((ITM_anti(i, 0) == 1) && (Payoff->operator()(Index_anti(i, t)) > C_hat_anti(i, 0)))
			{
				simulated_price_anti(i, 0) = Payoff->operator()(Index_anti(i, t));
				//early_exec.push_back(i);
			}

			{


				simulated_price_anti(i, 0) = df * simulated_price_anti(i, 0);
			}

			average_price(i, 0) = (simulated_price_anti(i,0) + simulated_price(i, 0)) * 0.5;
		}
	}

	MC_price = df*average_price.mean();

	MC_variance = average_price.variance();

};

AmericanMonteCarlo_basket_Antithetic_CV::AmericanMonteCarlo_basket_Antithetic_CV(size_t nbSimu, PayOffBasket* InputPayoff, PayOffBasket* CPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial, double closed_form_price):
	AmericanMonteCarlo(nbSimu, InputPayoff, diffusion, polynomial), x_diffusion(diffusion),CPayoff(CPayoff),ExpPriceClsForm(closed_form_price)
{
	r = diffusion->Get_rate();
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu / 2, 1);
	simulated_price_anti.Resize(nbSimu / 2, 1);
	average_price.Resize(nbSimu / 2, 1);
	simulated_price_CP.Resize(nbSimu / 2, 1);
	simulated_price_CP_anti.Resize(nbSimu / 2, 1);
};

void AmericanMonteCarlo_basket_Antithetic_CV::Simulate(double start, double end, size_t steps)
{

	std::cout << "American MC with antithetic and CV" << std::endl;

	double temp = 0.;
	double temp_anti = 0.;

	simulated_price.Resize(m_Simulation / 2, 1);
	simulated_price_anti.Resize(m_Simulation / 2, 1);
	average_price.Resize(m_Simulation / 2, 1);

	simulated_price_CP.Resize(m_Simulation / 2, 1);
	simulated_price_CP_anti.Resize(m_Simulation / 2, 1);
	average_price.Resize(m_Simulation / 2, 1);

	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation / 2, steps);
	matrix ITM(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix It(m_Simulation / 2, 1);
	matrix Phit(m_Simulation / 2, Phi.size());
	matrix poly_order_n(m_Simulation / 2, 1);
	matrix Phit_T(Phi.size(), m_Simulation / 2);
	matrix SDP_Phi(Phi.size(), Phi.size());
	matrix Phi_V(Phi.size(), 1);
	matrix inv_SDP_Phi(Phi.size(), Phi.size());
	matrix Beta_hat(Phi.size(), 1);
	matrix C_hat(m_Simulation / 2, 1);
	matrix chm(1, steps);
	matrix V(m_Simulation / 2, 1);
	matrix paths;
	matrix paths_anti;


	matrix chm_anti(1, steps);
	matrix Index_anti(m_Simulation / 2, steps);
	matrix ITM_anti(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix It_anti(m_Simulation / 2, 1);
	matrix Phit_anti(m_Simulation / 2, Phi.size());
	matrix poly_order_n_anti(m_Simulation / 2, 1);
	matrix Phit_T_anti(Phi.size(), m_Simulation / 2);
	matrix SDP_Phi_anti(Phi.size(), Phi.size());
	matrix Phi_V_anti(Phi.size(), 1);
	matrix inv_SDP_Phi_anti(Phi.size(), Phi.size());
	matrix Beta_hat_anti(Phi.size(), 1);
	matrix C_hat_anti(m_Simulation / 2, 1);
	matrix V_anti(m_Simulation / 2, 1);



	matrix log_Index(m_Simulation/2, steps);
	matrix ITM2(m_Simulation/2, 1); //need to separate which paths are ITM at each iteration 
	matrix It2(m_Simulation/2, 1);
	matrix Phit_2(m_Simulation/2, Phi.size());
	matrix poly_order_n_2(m_Simulation/2, 1);
	matrix Phit_T_2(Phi.size(), m_Simulation/2);
	matrix SDP_Phi_2(Phi.size(), Phi.size());
	matrix Phi_V_2(Phi.size(), 1);
	matrix inv_SDP_Phi_2(Phi.size(), Phi.size());
	matrix Beta_hat_2(Phi.size(), 1);
	matrix C_hat_2(m_Simulation/2, 1);
	matrix chm_2(1, steps);
	matrix V_2(m_Simulation/2, 1);
	matrix log_spot;


	matrix log_Index_anti(m_Simulation/2, steps);
	matrix ITM2_anti(m_Simulation/2, 1); //need to separate which paths are ITM at each iteration 
	matrix It2_anti(m_Simulation/2, 1);
	matrix Phit_2_anti(m_Simulation/2, Phi.size());
	matrix poly_order_n_2_anti(m_Simulation/2, 1);
	matrix Phit_T_2_anti(Phi.size(), m_Simulation/2);
	matrix SDP_Phi_2_anti(Phi.size(), Phi.size());
	matrix Phi_V_2_anti(Phi.size(), 1);
	matrix inv_SDP_Phi_2_anti(Phi.size(), Phi.size());
	matrix Beta_hat_2_anti(Phi.size(), 1);
	matrix C_hat_2_anti(m_Simulation/2, 1);
	matrix chm_2_anti(1, steps);
	matrix V_2_anti(m_Simulation/2, 1);
	matrix log_spot_anti;
	size_t k = 0;
	size_t k_a = 0;

	for (size_t s = 0; s < m_Simulation; s++)
	{

		//std::cout << "simulation " << s << std::endl;
		if (s % 2 == 0)
		{
			x_diffusion->Simulate(start, end, steps);
			paths = x_diffusion->GetAllPaths();
			log_spot = paths;
			chm = weights * paths;



			for (size_t i = 0; i < paths.nb_rows(); i++)
			{
				for (size_t j = 0; j < paths.nb_cols(); j++)
				{
					log_spot(i, j) = log(paths(i, j));
				}
			}

			chm_2 = weights * log_spot;

			for (size_t i = 0; i < steps; i++)
			{
				Index(k, i) = chm(0, i);
				log_Index(k, i) = chm_2(0, i);

			}

			simulated_price(k, 0) = Payoff->operator()(Index(k, Index.nb_cols() - 1));
			simulated_price_CP(k, 0) = CPayoff->operator()(log_Index(k, log_Index.nb_cols() - 1));

			k += 1;

		}
		else
		{
			paths_anti = x_diffusion->GetAllPathsAnti();
			log_spot_anti = paths_anti;
			chm_anti = weights * paths_anti;


			for (size_t i = 0; i < paths.nb_rows(); i++)
			{
				for (size_t j = 0; j < paths.nb_cols(); j++)
				{
					log_spot_anti(i, j) = log(paths_anti(i, j));
				}
			}

			chm_2_anti = weights * log_spot_anti;

			for (size_t i = 0; i < steps; i++)
			{
				Index_anti(k_a, i) = chm_anti(0, i);
				log_Index_anti(k_a, i) = chm_2_anti(0, i);

			}


			simulated_price_anti(k_a, 0) = Payoff->operator()(Index_anti(k_a, Index_anti.nb_cols() - 1));
			simulated_price_CP_anti(k_a,0) = CPayoff->operator()(log_Index_anti(k_a, log_Index_anti.nb_cols() - 1));
			k_a += 1;
		}
	}
	double df = exp(-r * x_diffusion->Get_Dt());

	// Regression part 

	for (size_t t = steps - 2; t > 0; t--)
	{


		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			V(i, 0) = df * simulated_price(i, 0);
			It(i, 0) = Index(i, t); //take index at t-1;

			V_2(i, 0) = df * simulated_price_CP(i, 0);
			It2(i, 0) = log_Index(i, t); //take index at t-1;

			V_anti(i, 0) = df * simulated_price_anti(i, 0);
			It_anti(i, 0) = Index_anti(i, t); //take index at t-1;

			V_2_anti(i, 0) = df * simulated_price_CP_anti(i, 0);
			It2_anti(i, 0) = log_Index_anti(i, t); //take index at t-1;


			if (Payoff->operator()(It(i, 0)) > 0)
			{
				ITM(i, 0) = 1;
			}
			else { ITM(i, 0) = -1; };

			if (CPayoff->operator()(It2(i, 0)) > 0)
			{
				ITM2(i, 0) = 1;
			}
			else { ITM2(i, 0) = -1; };


			if (Payoff->operator()(It_anti(i, 0)) > 0)
			{
				ITM_anti(i, 0) = 1;
			}
			else { ITM_anti(i, 0) = -1; };


			if (CPayoff->operator()(It2_anti(i, 0)) > 0)
			{
				ITM2_anti(i, 0) = 1;
			}
			else { ITM2_anti(i, 0) = -1; };

		}


		for (size_t b = 0; b < Phi.size(); b++)

		{
			poly_order_n = Phi[b]->operator()(It);
			poly_order_n_anti = Phi[b]->operator()(It_anti);

			poly_order_n_2 = Phi[b]->operator()(It2);
			poly_order_n_2_anti = Phi[b]->operator()(It2_anti);

			for (size_t i = 0; i < poly_order_n.nb_rows(); i++)
			{
				Phit(i, b) = poly_order_n(i, 0);

				Phit_anti(i, b) = poly_order_n_anti(i, 0);


				Phit_2(i, b) = poly_order_n_2(i, 0);

				Phit_2_anti(i, b) = poly_order_n_2_anti(i, 0);
			}


		}

		Phit_T = transpose(Phit);
		Phit_T_anti = transpose(Phit_anti);
		Phit_T_2 = transpose(Phit_2);
		Phit_T_2_anti = transpose(Phit_2_anti);

		SDP_Phi = Phit_T * Phit;
		SDP_Phi_anti = Phit_T_anti * Phit_anti;
		SDP_Phi_2 = Phit_T_2 * Phit_2;
		SDP_Phi_2_anti = Phit_T_2_anti * Phit_2_anti;


		Phi_V = Phit_T * V;
		Phi_V_anti = Phit_T_anti * V_anti;
		Phi_V_2 = Phit_T_2 * V_2;
		Phi_V_2_anti = Phit_T_2_anti * V_2_anti;


		inv_SDP_Phi = Inverse(SDP_Phi, SDP_Phi.nb_cols());
		inv_SDP_Phi_anti = Inverse(SDP_Phi_anti, SDP_Phi_anti.nb_cols());
		inv_SDP_Phi_2 = Inverse(SDP_Phi_2, SDP_Phi_2.nb_cols());
		inv_SDP_Phi_2_anti = Inverse(SDP_Phi_2_anti, SDP_Phi_2_anti.nb_cols());


		Beta_hat = inv_SDP_Phi * Phi_V;
		Beta_hat_anti = inv_SDP_Phi_anti * Phi_V_anti;
		Beta_hat_2 = inv_SDP_Phi_2 * Phi_V_2;
		Beta_hat_2_anti = inv_SDP_Phi_2_anti * Phi_V_2_anti;

		C_hat = Phit * Beta_hat;
		C_hat_anti = Phit_anti * Beta_hat_anti;
		C_hat_2 = Phit_2 * Beta_hat_2;
		C_hat_2_anti = Phit_2_anti * Beta_hat_2_anti;

		for (size_t i = 0; i < C_hat.nb_rows(); i++)
		{

			if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(i, t)) > C_hat(i, 0)))
			{
				simulated_price(i, 0) = Payoff->operator()(Index(i, t));
				//early_exec.push_back(i);
			}
			else
			{

				simulated_price(i, 0) = df * simulated_price(i, 0);
			}

			if ((ITM2(i, 0) == 1) && (CPayoff->operator()(log_Index(i, t)) > C_hat_2(i, 0)))
			{
				simulated_price_CP(i, 0) = CPayoff->operator()(log_Index(i, t));
				//early_exec.push_back(i);
			}
			else
			{

				simulated_price_CP(i, 0) = df * simulated_price_CP(i, 0);
			}


			if ((ITM_anti(i, 0) == 1) && (Payoff->operator()(Index_anti(i, t)) > C_hat_anti(i, 0)))
			{
				simulated_price_anti(i, 0) = Payoff->operator()(Index_anti(i, t));
				//early_exec.push_back(i);
			}
			else
			{

				simulated_price_anti(i, 0) = df * simulated_price_anti(i, 0);
			}


			if ((ITM2_anti(i, 0) == 1) && (CPayoff->operator()(log_Index_anti(i, t)) > C_hat_2_anti(i, 0)))
			{
				simulated_price_CP_anti(i, 0) = CPayoff->operator()(log_Index_anti(i, t));
				//early_exec.push_back(i);
			}
			else
			{

				simulated_price_CP_anti(i, 0) = df * simulated_price_CP_anti(i, 0);
			}

			temp = simulated_price(i, 0) - simulated_price_CP(i,0);
			temp_anti = simulated_price_anti(i,0) - simulated_price_CP_anti(i,0);

			average_price(i, 0) = (temp + temp_anti) * 0.5 + ExpPriceClsForm;

			temp = 0.;
			temp_anti = 0;
		}
	}

	MC_price = df*average_price.mean();

	MC_variance = average_price.variance();

};

/////////////////// BERMUDEAN ///////////////////////////////////////////////////////

BermudeanMonteCarlo::BermudeanMonteCarlo(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial, CalendarManagement* wkday, matrix exec_schedule):
	AmericanMonteCarlo(nbSimu,InputPayoff,diffusion,polynomial), 
	wkday(wkday), exec_schedule(exec_schedule)
{	
	r = diffusion->Get_rate();
	MC_price = 0.;
	MC_variance = 0.;
	//construction of the Dt_Schedule matrix that will be used to interpolate the spot in execution date
	Dt_schedule.Resize(exec_schedule.nb_rows() - 1, 1);


};


Bermudean_BasketOption::Bermudean_BasketOption(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial, CalendarManagement* wkday, matrix exec_schedule):
	BermudeanMonteCarlo(nbSimu, InputPayoff, diffusion,polynomial,wkday,exec_schedule)
{
	r = diffusion->Get_rate();
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu, 1);
};



void Bermudean_BasketOption::Simulate(double start, double end, size_t steps)
{
		//matrix Phit(m_Simulation, Phi.size());
		//matrix poly_order_n(m_Simulation, 1);
		//matrix Phit_T(Phi.size(), m_Simulation);
		//matrix SDP_Phi(Phi.size(), Phi.size());
		//matrix Phi_V(Phi.size(), 1);
		//matrix inv_SDP_Phi(Phi.size(), Phi.size());
		//matrix Beta_hat(Phi.size(), 1);


	std::cout << "Bermudean BasketOption LS" << std::endl;
	simulated_price.Resize(m_Simulation, 1);
	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation, exec_schedule.nb_rows());
	matrix ITM(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix It(m_Simulation, 1);
	matrix C_hat(m_Simulation, 1);
	matrix chm(1, exec_schedule.nb_rows());
	matrix V(m_Simulation, 1);
	matrix paths;
	matrix interpolated_paths;
	double dt_sde;
	double df;
	//double test_df;

	for (size_t s = 0; s < m_Simulation; s++)
	{
		//std::cout << "simulation " << s << std::endl;
		m_diffusion->Simulate(start, end, steps);
		paths = m_diffusion->GetAllPaths();
		interpolated_paths = wkday->index_executed(paths, exec_schedule, steps);

		

		//interpolated_paths.Print();
		chm = weights * interpolated_paths;
		//std::cout << "end simulation " << s << std::endl;
		for (size_t i = 0; i < chm.nb_cols(); i++)
		{
			Index(s, i) = chm(0, i);

		}
		//std::cout << "first value of chm  at simulation " <<  chm(0,0) << std::endl;
		//std::cout << "4th value of chm  at simulation " << chm(0, 4) << std::endl;
		//std::cout << "last value of chm  at simulation " << chm(0,9) << std::endl;


		simulated_price(s, 0) = Payoff->operator()(Index(s, Index.nb_cols() - 1));
		//populate the vector at maturity of the payoff

	}

	dt_sde = m_diffusion->Get_Dt();

	for (size_t j = 1; j < exec_schedule.nb_rows(); j++)
	{

		Dt_schedule(j - 1, 0) = (ceil(exec_schedule(j, 0) * steps) - ceil(exec_schedule(j - 1, 0) * steps)) * dt_sde;

	}

	double last_dt = ceil(exec_schedule(0, 0) * steps) * dt_sde;
	// Regression part 

	//std::cout << " start regression" << std::endl;


	size_t Isteps = Index.nb_cols() -2; //Index will have as much columns as there are execution dates.


	//Dt_schedule.Print();
	//std::cout << Isteps << std::endl;
	for (size_t t = Isteps; t > 0; t--)
	{

		df = exp(-r * Dt_schedule(t, 0));
		//df = exp(-r * Dt_schedule(t, 0));
		//test_df = exp(-r * dt_sde);
		//std::cout << "constant dt " << test_df << std::endl;
		//std::cout << " dt schedule " << df << std::endl;
		//std::cout << " Dt " << Dt_schedule(t, 0) << std::endl;
		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			V(i, 0) = df * simulated_price(i, 0);
			It(i, 0) = Index(i, t); //take index at t-1;

			if (Payoff->operator()(It(i, 0)) > 0)
			{
				ITM(i, 0) = 1;
			}
			else { ITM(i, 0) = -1; };
		}


		C_hat = C_Hat_regression(It, V); //the function will apply the polynomes, create Beta_Hat, and return C_hat as a matrix

		for (size_t i = 0; i < C_hat.nb_rows(); i++)
		{

			if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(i, t)) > C_hat(i, 0)))
			{
				simulated_price(i, 0) = Payoff->operator()(Index(i, t));
				//early_exec.push_back(i);
			}

			{
				simulated_price(i, 0) = df * simulated_price(i, 0);
			}


		}

	}

	
	double last_df = exp(-r *last_dt);
	MC_price = last_df*simulated_price.mean();

	MC_variance = simulated_price.variance();

};

Bermudean_BasketOption_CV::Bermudean_BasketOption_CV(size_t nbSimu, PayOffBasket* InputPayoff, PayOffBasket* control_payoff,RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial, CalendarManagement* wkday, matrix exec_schedule, double closed_form_price)
	:BermudeanMonteCarlo(nbSimu, InputPayoff, diffusion, polynomial, wkday, exec_schedule),
	CPayoff(control_payoff), ExpPriceClsForm(closed_form_price)
{
	r = diffusion->Get_rate();
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu, 1);
	simulated_price_CP.Resize(nbSimu, 1); //to compute the mean of the control variate payoff 
};


void Bermudean_BasketOption_CV::Simulate(double start, double end, size_t steps)
{
	//matrix Phit(m_Simulation, Phi.size());
	//matrix poly_order_n(m_Simulation, 1);
	//matrix Phit_T(Phi.size(), m_Simulation);
	//matrix SDP_Phi(Phi.size(), Phi.size());
	//matrix Phi_V(Phi.size(), 1);
	//matrix inv_SDP_Phi(Phi.size(), Phi.size());
	//matrix Beta_hat(Phi.size(), 1);

	//matrix Phit_2(m_Simulation, Phi.size());
	//matrix poly_order_n_2(m_Simulation, 1);
	//matrix Phit_T_2(Phi.size(), m_Simulation);
	//matrix SDP_Phi_2(Phi.size(), Phi.size());
	//matrix Phi_V_2(Phi.size(), 1);
	//matrix inv_SDP_Phi_2(Phi.size(), Phi.size());
	//matrix Beta_hat_2(Phi.size(), 1);

	std::cout << "Bermudean MC LS with Control Variable" << std::endl;

	simulated_price.Resize(m_Simulation, 1);
	simulated_price_CP.Resize(m_Simulation, 1);
	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation, exec_schedule.nb_rows());
	matrix ITM(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix It(m_Simulation, 1);
	matrix C_hat(m_Simulation, 1);
	matrix chm;
	matrix V(m_Simulation, 1);
	matrix paths;

	matrix log_Index(m_Simulation, exec_schedule.nb_rows());
	matrix ITM2(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix It2(m_Simulation, 1);
	matrix C_hat_2(m_Simulation, 1);
	matrix chm_2;
	matrix V_2(m_Simulation, 1);
	matrix log_spot;
	matrix interpolated_paths;
	double dt_sde;
	double df;

	matrix CV(m_Simulation, 1);

	for (size_t s = 0; s < m_Simulation; s++)
	{
		m_diffusion->Simulate(start, end, steps);
		paths = m_diffusion->GetAllPaths();
		log_spot = paths;
		interpolated_paths = wkday->index_executed(paths, exec_schedule,steps);

		for (size_t i = 0; i < interpolated_paths.nb_rows(); i++)
		{
			for (size_t j = 0; j < interpolated_paths.nb_cols(); j++)
			{
				log_spot(i, j) = log(interpolated_paths(i, j));
			}
		}
		chm = weights * interpolated_paths;
		chm_2 = weights * log_spot;

		for (size_t i = 0; i < interpolated_paths.nb_cols(); i++)
		{
			Index(s, i) = chm(0, i);
			log_Index(s, i) = chm_2(0, i);
		}

		simulated_price(s, 0) = Payoff->operator()(Index(s, Index.nb_cols() - 1));
		simulated_price_CP(s, 0) = CPayoff->operator()(log_Index(s, Index.nb_cols() - 1));
		//populate the vector at maturity of the payoff

	}

	dt_sde = m_diffusion->Get_Dt();

	for (size_t j = 1; j < exec_schedule.nb_rows(); j++)
	{

		Dt_schedule(j - 1, 0) = (ceil(exec_schedule(j, 0) * steps) - ceil(exec_schedule(j - 1, 0) * steps)) * dt_sde;

	}

	double last_dt = ceil(exec_schedule(0, 0) * steps) * dt_sde;

	// Regression part

	size_t Isteps = Index.nb_cols() - 2;

	for (size_t t = Isteps; t > 0; t--)
	{
		df = exp(-r * Dt_schedule(t, 0));
		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			V(i, 0) = df * simulated_price(i, 0);
			V_2(i, 0) = df * simulated_price_CP(i, 0);
			It(i, 0) = Index(i, t); //take index at t-1;
			It2(i, 0) = log_Index(i, t);
			if (Payoff->operator()(It(i, 0)) > 0)
			{
				ITM(i, 0) = 1;
			}
			else { ITM(i, 0) = -1; };

			if (CPayoff->operator()(It2(i, 0)) > 0)
			{
				ITM2(i, 0) = 1;
			}
			else { ITM2(i, 0) = -1; };
		}

		C_hat = C_Hat_regression(It, V);
		C_hat_2 = C_Hat_regression(It2, V_2);

		for (size_t i = 0; i < C_hat.nb_rows(); i++)
		{

			if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(i, t)) > C_hat(i, 0)))
			{
				simulated_price(i, 0) = Payoff->operator()(Index(i, t));
			}

			{
				simulated_price(i, 0) = df * simulated_price(i, 0);
			}




			if ((ITM2(i, 0) == 1) && (CPayoff->operator()(log_Index(i, t)) > C_hat_2(i, 0)))
			{
				simulated_price_CP(i, 0) = CPayoff->operator()(log_Index(i, t));
			}

			{


				simulated_price_CP(i, 0) = df * simulated_price_CP(i, 0);
			}

			CV(i, 0) = simulated_price(i, 0) - simulated_price_CP(i, 0) + ExpPriceClsForm;;


		}
	}

	double last_df = exp(-r * last_dt);
	MC_price = last_df*CV.mean();
	MC_variance = CV.variance();

};

Bermudean_BasketOption_antithetic::Bermudean_BasketOption_antithetic(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial, CalendarManagement* wkday, matrix exec_schedule)
	:BermudeanMonteCarlo(nbSimu, InputPayoff, diffusion, polynomial, wkday, exec_schedule), x_diffusion(diffusion)
{
	r = diffusion->Get_rate();
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu / 2, 1);
	simulated_price_anti.Resize(nbSimu / 2, 1);
	average_price.Resize(nbSimu / 2, 1);

};

void Bermudean_BasketOption_antithetic::Simulate(double start, double end, size_t steps)
{
	//matrix Phit(m_Simulation / 2, Phi.size());
	//matrix poly_order_n(m_Simulation / 2, 1);
	//matrix Phit_T(Phi.size(), m_Simulation / 2);
	//matrix SDP_Phi(Phi.size(), Phi.size());
	//matrix Phi_V(Phi.size(), 1);
	//matrix inv_SDP_Phi(Phi.size(), Phi.size());
	//matrix Beta_hat(Phi.size(), 1);

	//matrix Phit(m_Simulation / 2, Phi.size());
	//matrix poly_order_n(m_Simulation / 2, 1);
	//matrix Phit_T(Phi.size(), m_Simulation / 2);
	//matrix SDP_Phi(Phi.size(), Phi.size());
	//matrix Phi_V(Phi.size(), 1);
	//matrix inv_SDP_Phi(Phi.size(), Phi.size());
	//matrix Beta_hat(Phi.size(), 1);

	std::cout << "Bermudean MC with antithetic" << std::endl;

	simulated_price.Resize(m_Simulation / 2, 1);
	simulated_price_anti.Resize(m_Simulation / 2, 1);
	average_price.Resize(m_Simulation / 2, 1);

	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation / 2, exec_schedule.nb_rows());
	matrix ITM(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix It(m_Simulation / 2, 1);
	matrix C_hat(m_Simulation / 2, 1);
	matrix chm;
	matrix V(m_Simulation / 2, 1);
	matrix paths;
	matrix interpolated_paths;
	matrix paths_anti;
	matrix interpolated_paths_anti;
	

	matrix chm_anti;
	matrix Index_anti(m_Simulation / 2, exec_schedule.nb_rows());
	matrix ITM_anti(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix It_anti(m_Simulation / 2, 1);
	matrix C_hat_anti(m_Simulation / 2, 1);
	matrix V_anti(m_Simulation / 2, 1);
	double dt_sde;
	double df;
	//double test_df;
	size_t k = 0;
	size_t k_a = 0;

	for (size_t s = 0; s < m_Simulation; s++)
	{

		//std::cout << "simulation " << s << std::endl;
		if (s % 2 == 0)
		{
			x_diffusion->Simulate(start, end, steps);
			paths = x_diffusion->GetAllPaths();
			interpolated_paths = wkday->index_executed(paths, exec_schedule, steps);
			chm = weights * interpolated_paths;

			for (size_t i = 0; i < chm.nb_cols(); i++)
			{
				Index(k, i) = chm(0, i);

			}

			simulated_price(k, 0) = Payoff->operator()(Index(k, Index.nb_cols() - 1));

			k += 1;



		}
		else
		{
			paths_anti = x_diffusion->GetAllPathsAnti();
			interpolated_paths_anti = wkday->index_executed(paths_anti, exec_schedule, steps);
			chm_anti = weights * interpolated_paths_anti;

			for (size_t i = 0; i < chm_anti.nb_cols(); i++)
			{
				Index_anti(k_a, i) = chm_anti(0, i);

			}

			simulated_price_anti(k_a, 0) = Payoff->operator()(Index_anti(k_a, Index_anti.nb_cols() - 1));
			k_a += 1;
		}
	}

	dt_sde = x_diffusion->Get_Dt();

	for (size_t j = 1; j < exec_schedule.nb_rows(); j++)
	{

		Dt_schedule(j - 1, 0) = (ceil(exec_schedule(j, 0) * steps) - ceil(exec_schedule(j - 1, 0) * steps)) * dt_sde;

	}

	double last_dt = ceil(exec_schedule(0, 0) * steps) * dt_sde;


	// Regression part 
	size_t Isteps = Index.nb_cols() - 2;

	for (size_t t = Isteps; t > 0; t--)
	{

		df = exp(-r * Dt_schedule(t, 0));
		//test_df = exp(-r * dt_sde);
		//std::cout << "constant dt " << test_df << std::endl;
		//std::cout << " dt schedule " << df << std::endl;
		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			V(i, 0) = df * simulated_price(i, 0);
			It(i, 0) = Index(i, t); //take index at t-1;

			if (Payoff->operator()(It(i, 0)) > 0)
			{
				ITM(i, 0) = 1;
			}
			else { ITM(i, 0) = -1; };

			V_anti(i, 0) = df * simulated_price_anti(i, 0);
			It_anti(i, 0) = Index_anti(i, t); //take index at t-1;

			if (Payoff->operator()(It_anti(i, 0)) > 0)
			{
				ITM_anti(i, 0) = 1;
			}
			else { ITM_anti(i, 0) = -1; };
		}

		
		C_hat = C_Hat_regression(It, V);

		C_hat_anti = C_Hat_regression(It_anti, V_anti);


		for (size_t i = 0; i < C_hat.nb_rows(); i++)
		{

			if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(i, t)) > C_hat(i, 0)))
			{
				simulated_price(i, 0) = Payoff->operator()(Index(i, t));
				//early_exec.push_back(i);
			}

			{
				simulated_price(i, 0) = df * simulated_price(i, 0);
			}


			if ((ITM_anti(i, 0) == 1) && (Payoff->operator()(Index_anti(i, t)) > C_hat_anti(i, 0)))
			{
				simulated_price_anti(i, 0) = Payoff->operator()(Index_anti(i, t));
				//early_exec.push_back(i);
			}

			{

				simulated_price_anti(i, 0) = df * simulated_price_anti(i, 0);
			}

			average_price(i, 0) = (simulated_price_anti(i, 0) + simulated_price(i, 0)) * 0.5;
		}
	}

	double last_df = exp(-r * last_dt);
	MC_price = last_df*average_price.mean();

	MC_variance = average_price.variance();

};

Bermudean_BasketOption_antithetic_CV::Bermudean_BasketOption_antithetic_CV(size_t nbSimu, PayOffBasket* InputPayoff, PayOffBasket* control_payoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial,CalendarManagement* wkday, matrix exec_schedule, double closed_form_price):
	BermudeanMonteCarlo(nbSimu, InputPayoff, diffusion, polynomial, wkday, exec_schedule),
	x_diffusion(diffusion), CPayoff(control_payoff), ExpPriceClsForm(closed_form_price)
{
	r = diffusion->Get_rate();
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu / 2, 1);
	simulated_price_anti.Resize(nbSimu / 2, 1);
	average_price.Resize(nbSimu / 2, 1);

	simulated_price_CP.Resize(nbSimu / 2, 1);
	simulated_price_CP_anti.Resize(nbSimu / 2, 1);
	average_price.Resize(nbSimu / 2, 1);
};

void Bermudean_BasketOption_antithetic_CV::Simulate(double start, double end, size_t steps)
{

	std::cout << "Bermudean MC with antithetic and CV" << std::endl;

	double temp = 0.;
	double temp_anti = 0.;

	simulated_price.Resize(m_Simulation / 2, 1);
	simulated_price_anti.Resize(m_Simulation / 2, 1);
	average_price.Resize(m_Simulation / 2, 1);

	simulated_price_CP.Resize(m_Simulation / 2, 1);
	simulated_price_CP_anti.Resize(m_Simulation / 2, 1);
	average_price.Resize(m_Simulation / 2, 1);

	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation / 2, exec_schedule.nb_rows());
	matrix ITM(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix It(m_Simulation / 2, 1);
	//matrix Phit(m_Simulation / 2, Phi.size());
	//matrix poly_order_n(m_Simulation / 2, 1);
	//matrix Phit_T(Phi.size(), m_Simulation / 2);
	//matrix SDP_Phi(Phi.size(), Phi.size());
	//matrix Phi_V(Phi.size(), 1);
	//matrix inv_SDP_Phi(Phi.size(), Phi.size());
	//matrix Beta_hat(Phi.size(), 1);
	matrix C_hat(m_Simulation / 2, 1);
	matrix chm;
	matrix V(m_Simulation / 2, 1);
	matrix paths;
	matrix paths_anti;


	matrix chm_anti;
	matrix Index_anti(m_Simulation / 2, exec_schedule.nb_rows());
	matrix ITM_anti(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix It_anti(m_Simulation / 2, 1);
	//matrix Phit_anti(m_Simulation / 2, Phi.size());
	//matrix poly_order_n_anti(m_Simulation / 2, 1);
	//matrix Phit_T_anti(Phi.size(), m_Simulation / 2);
	//matrix SDP_Phi_anti(Phi.size(), Phi.size());
	//matrix Phi_V_anti(Phi.size(), 1);
	//matrix inv_SDP_Phi_anti(Phi.size(), Phi.size());
	//matrix Beta_hat_anti(Phi.size(), 1);
	matrix C_hat_anti(m_Simulation / 2, 1);
	matrix V_anti(m_Simulation / 2, 1);



	matrix log_Index(m_Simulation, exec_schedule.nb_rows());
	matrix ITM2(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix It2(m_Simulation, 1);
	//matrix Phit_2(m_Simulation, Phi.size());
	//matrix poly_order_n_2(m_Simulation, 1);
	//matrix Phit_T_2(Phi.size(), m_Simulation);
	//matrix SDP_Phi_2(Phi.size(), Phi.size());
	//matrix Phi_V_2(Phi.size(), 1);
	//matrix inv_SDP_Phi_2(Phi.size(), Phi.size());
	//matrix Beta_hat_2(Phi.size(), 1);
	matrix C_hat_2(m_Simulation, 1);
	matrix chm_2;
	matrix V_2(m_Simulation, 1);
	matrix log_spot;


	matrix log_Index_anti(m_Simulation, exec_schedule.nb_rows());
	matrix ITM2_anti(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix It2_anti(m_Simulation, 1);
	//matrix Phit_2_anti(m_Simulation, Phi.size());
	//matrix poly_order_n_2_anti(m_Simulation, 1);
	//matrix Phit_T_2_anti(Phi.size(), m_Simulation);
	//matrix SDP_Phi_2_anti(Phi.size(), Phi.size());
	//matrix Phi_V_2_anti(Phi.size(), 1);
	//matrix inv_SDP_Phi_2_anti(Phi.size(), Phi.size());
	//matrix Beta_hat_2_anti(Phi.size(), 1);
	matrix C_hat_2_anti(m_Simulation, 1);
	matrix chm_2_anti;
	matrix V_2_anti(m_Simulation, 1);
	matrix log_spot_anti;

	matrix interpolated_paths;
	matrix interpolated_paths_anti;
	double dt_sde;
	double df;

	size_t k = 0;
	size_t k_a = 0;

	for (size_t s = 0; s < m_Simulation; s++)
	{

		std::cout << "simulation " << s << std::endl;
		if (s % 2 == 0)
		{
			x_diffusion->Simulate(start, end, steps);
			paths = x_diffusion->GetAllPaths();
			interpolated_paths = wkday->index_executed(paths, exec_schedule, steps);
			log_spot = paths;
			chm = weights * interpolated_paths;

			for (size_t i = 0; i < interpolated_paths.nb_rows(); i++)
			{
				for (size_t j = 0; j < interpolated_paths.nb_cols(); j++)
				{
					log_spot(i, j) = log(interpolated_paths(i, j));
				}
			}

			chm_2 = weights * log_spot;

			for (size_t i = 0; i < chm.nb_cols(); i++)
			{
				Index(k, i) = chm(0, i);
				log_Index(k, i) = chm_2(0, i);

			}


			simulated_price(k, 0) = Payoff->operator()(Index(k, Index.nb_cols() - 1));
			simulated_price_CP(k, 0) = CPayoff->operator()(log_Index(k, log_Index.nb_cols() - 1));

			k += 1;
		}
		else
		{
			paths_anti = x_diffusion->GetAllPathsAnti();
			interpolated_paths_anti = wkday->index_executed(paths_anti, exec_schedule, steps);
			log_spot_anti = interpolated_paths_anti;
			chm_anti = weights * interpolated_paths_anti;


			for (size_t i = 0; i < interpolated_paths_anti.nb_rows(); i++)
			{
				for (size_t j = 0; j < interpolated_paths_anti.nb_cols(); j++)
				{
					log_spot_anti(i, j) = log(interpolated_paths_anti(i, j));
				}
			}

			chm_2_anti = weights * log_spot_anti;

			for (size_t i = 0; i < chm_anti.nb_cols(); i++)
			{
				Index_anti(k_a, i) = chm_anti(0, i);
				log_Index_anti(k_a, i) = chm_2_anti(0, i);

			}


			simulated_price_anti(k_a, 0) = Payoff->operator()(Index_anti(k_a, Index_anti.nb_cols() - 1));
			simulated_price_CP_anti(k_a, 0) = Payoff->operator()(log_Index_anti(k_a, log_Index_anti.nb_cols() - 1));
			k_a += 1;
		}
	}

	dt_sde = x_diffusion->Get_Dt();

	for (size_t j = 1; j < exec_schedule.nb_rows(); j++)
	{

		Dt_schedule(j - 1, 0) = (ceil(exec_schedule(j, 0)*steps) - ceil(exec_schedule(j - 1, 0)*steps)) * dt_sde;

	}

	double last_dt = ceil(exec_schedule(0, 0)) * dt_sde;

	// Regression part 

	size_t Isteps = Index.nb_cols() - 2;

	for (size_t t = Isteps; t > 0; t--)
	{

		df = exp(-r * Dt_schedule(t, 0));
		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			V(i, 0) = df * simulated_price(i, 0);
			It(i, 0) = Index(i, t); //take index at t-1;

			V_2(i, 0) = df * simulated_price_CP(i, 0);
			It2(i, 0) = log_Index(i, t); //take index at t-1;

			V_anti(i, 0) = df * simulated_price_anti(i, 0);
			It_anti(i, 0) = Index_anti(i, t); //take index at t-1;

			V_2_anti(i, 0) = df * simulated_price_CP_anti(i, 0);
			It2_anti(i, 0) = log_Index_anti(i, t); //take index at t-1;


			if (Payoff->operator()(It(i, 0)) > 0)
			{
				ITM(i, 0) = 1;
			}
			else { ITM(i, 0) = -1; };

			if (CPayoff->operator()(It2(i, 0)) > 0)
			{
				ITM2(i, 0) = 1;
			}
			else { ITM2(i, 0) = -1; };


			if (Payoff->operator()(It_anti(i, 0)) > 0)
			{
				ITM_anti(i, 0) = 1;
			}
			else { ITM_anti(i, 0) = -1; };


			if (CPayoff->operator()(It2_anti(i, 0)) > 0)
			{
				ITM2_anti(i, 0) = 1;
			}
			else { ITM2_anti(i, 0) = -1; };

		}


		C_hat = C_Hat_regression(It, V);
		C_hat_2 = C_Hat_regression(It2, V_2);
		C_hat_anti = C_Hat_regression(It_anti, V_anti);
		C_hat_2_anti = C_Hat_regression(It2_anti, V_2_anti);

		for (size_t i = 0; i < C_hat.nb_rows(); i++)
		{

			if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(i, t)) > C_hat(i, 0)))
			{
				simulated_price(i, 0) = Payoff->operator()(Index(i, t));
				//early_exec.push_back(i);
			}
			else
			{

				simulated_price(i, 0) = df * simulated_price(i, 0);
			}

			if ((ITM2(i, 0) == 1) && (CPayoff->operator()(log_Index(i, t)) > C_hat_2(i, 0)))
			{
				simulated_price_CP(i, 0) = CPayoff->operator()(log_Index(i, t));
				//early_exec.push_back(i);
			}
			else
			{

				simulated_price_CP(i, 0) = df * simulated_price_CP(i, 0);
			}


			if ((ITM_anti(i, 0) == 1) && (Payoff->operator()(Index_anti(i, t)) > C_hat_anti(i, 0)))
			{
				simulated_price_anti(i, 0) = Payoff->operator()(Index_anti(i, t));
				//early_exec.push_back(i);
			}
			else
			{

				simulated_price_anti(i, 0) = df * simulated_price_anti(i, 0);
			}


			if ((ITM2_anti(i, 0) == 1) && (CPayoff->operator()(log_Index_anti(i, t)) > C_hat_2_anti(i, 0)))
			{
				simulated_price_CP_anti(i, 0) = CPayoff->operator()(log_Index_anti(i, t));
				//early_exec.push_back(i);
			}
			else
			{

				simulated_price_CP_anti(i, 0) = df * simulated_price_CP_anti(i, 0);
			}

			temp = simulated_price(i, 0) - simulated_price_CP(i, 0) + ExpPriceClsForm;
			temp_anti = simulated_price_anti(i, 0) - simulated_price_CP_anti(i, 0) + ExpPriceClsForm;

			average_price(i, 0) = (temp + temp_anti) * 0.5;
		}
	}

	double last_df = exp(-r * last_dt);
	MC_price = last_df*average_price.mean();

	MC_variance = average_price.variance();

};