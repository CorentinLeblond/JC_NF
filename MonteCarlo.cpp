#include "MonteCarlo.h"
#include <iostream>

//MC
//////////////////////////////////////////////////////////

double MonteCarlo::GetPrice(double& r, double& T)
{
	return exp(-r * T) * MC_price;
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

	MC_price = simulated_price.mean();
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
	MC_price = simulated_price.mean();
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
<<<<<<< HEAD
	std::cout << "MC European Basket and Anti"<< std::endl;
	
	simulated_price.Resize(m_Simulation/2,1);
	// simulated_price_Anti.Resize(m_Simulation/2,1);
	size_t k = 0;
	
=======
	std::cout << "MC European Basket antithetic"<< std::endl;

>>>>>>> LSM working in vanilla, CV and antithetic
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
		}		
		k+=1;
	}

	MC_price = average_price.mean();
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
			simulated_price(k,0) = Payoff->operator()(maturity_spot)-CPayoff->operator()(maturity_spot)+ExpPriceClsForm;;
		}
		else
		{
			simulated_price_Anti(k,0) = Payoff->operator()(maturity_spot)- CPayoff->operator()(maturity_spot)+ExpPriceClsForm;;
			average_price(k,0) = 0.5*(simulated_price_Anti(k,0) + simulated_price(k,0));
		}		
		k+=1;
	}

<<<<<<< HEAD
	MC_price = average_price.mean();
	MC_variance = average_price.variance();
	
};	
=======
};

AmericanMonteCarlo::AmericanMonteCarlo(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial, double df) :
	Phi(polynomial), df(df)
{
	Payoff = InputPayoff;
	m_Simulation = nbSimu;
	m_diffusion = diffusion;
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu,1);
};


void AmericanMonteCarlo::Simulate(double start, double end, size_t steps)
{

	std::cout << "American MC LS" << std::endl;
	matrix weights = Payoff->GetWeights();
	//weights.Print();
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
		//paths.Print();
		chm = weights * paths;
		//chm.Print();
		//Index.addrows(chm.GetMatrix());
		//appendrow(Index, chm);//void function that modifies the Index matrix 
		//Index.Print();

		for (size_t i = 0; i < steps; i++) 
		{
			Index(s, i) = chm(0,i);
		
		}

		simulated_price(s, 0) = Payoff->operator()(Index(s,Index.nb_cols()-1));
		//simulated_price.Print();
		//populate the vector at maturity of the payoff

		//std::cout << "here" << std::endl;
	}

	// Regression part 

	//Index.Print();

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

		//2) construct Phi that contains the polynominial transformation at each order 

		//for (size_t i = 0; i < Index.nb_rows(); i++) 
		//{
		//	
		//	
		//			
		//}
		//It = Index.area(Index.nb_rows() - 1, t, { 0,t });
		//std::cout << "It" << std::endl;
		//It.Print();
		for (size_t b = 0; b < Phi.size(); b++) 
		
		{	
			poly_order_n = Phi[b]->operator()(It);
			//std::cout << "plynome " << b << std::endl;
			//poly_order_n.Print();

			for (size_t i = 0; i < poly_order_n.nb_rows(); i++)
			{
				Phit(i, b) = poly_order_n(i,0);
			}
			

		}

		//std::cout << "Phit " << std::endl;
		//Phit.Print();

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

	MC_price = simulated_price.mean();

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


AmericanMonteCarlo_controlevariable::AmericanMonteCarlo_controlevariable(size_t nbSimu, PayOffBasket* InputPayoff, PayOffBasket* CPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial, double df)
	:AmericanMonteCarlo(nbSimu,InputPayoff,diffusion,polynomial,df), CPayoff(CPayoff)
{
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu, 1);
	simulated_price_CP.Resize(nbSimu,1); //to compute the mean of the control variate payoff 
};


void AmericanMonteCarlo_controlevariable::Simulate(double start, double end, size_t steps)
{

	std::cout << "American MC LS with Control Variable" << std::endl;
	matrix weights = Payoff->GetWeights();
	//weights.Print();
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
		//std::cout << "simulation " << s << std::endl;
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
		//simulated_price.Print();
		//populate the vector at maturity of the payoff

		//std::cout << "here" << std::endl;
	}

	// Regression part 

	//Index.Print();

	for (size_t t = steps - 2; t > 0; t--)
	{

		//std::cout << "regression " << t << std::endl;
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

		//2) construct Phi that contains the polynominial transformation at each order 

		for (size_t b = 0; b < Phi.size(); b++)

		{
			poly_order_n = Phi[b]->operator()(It);
			poly_order_n_2 = Phi[b]->operator()(It2);
			//std::cout << "plynome " << b << std::endl;
			//poly_order_n.Print();

			for (size_t i = 0; i < poly_order_n.nb_rows(); i++)
			{
				Phit(i, b) = poly_order_n(i, 0);
				Phit_2(i, b) = poly_order_n_2(i, 0);
			}


		}

		//std::cout << "Phit " << std::endl;
		//Phit.Print();

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
				early_exec.push_back(i);
			}

			{


				simulated_price(i, 0) = df * simulated_price(i, 0);
			}




			if ((ITM2(i, 0) == 1) && (CPayoff->operator()(log_Index(i, t)) > C_hat_2(i, 0)))
			{
				simulated_price_CP(i, 0) = CPayoff->operator()(log_Index(i, t));
				//early_exec.push_back(i);
			}

			{


				simulated_price_CP(i, 0) = df * simulated_price_CP(i, 0);
			}

			CV(i, 0) = simulated_price(i, 0) - simulated_price_CP(i, 0);


		}
	} 


	MC_price = CV.mean() + simulated_price_CP.mean() ;

	MC_variance = CV.variance();

};

AmericanMonteCarlo_Antithetic::AmericanMonteCarlo_Antithetic(size_t nbSimu, PayOffBasket* InputPayoff, RandomProcess* diffusion,
	std::vector<basis_functions*> polynomial, double df) 
	:AmericanMonteCarlo(nbSimu,InputPayoff,diffusion,polynomial,df), x_diffusion(diffusion)
{
	MC_price = 0.;
	MC_variance = 0.;
	simulated_price.Resize(nbSimu, 1);

};

void AmericanMonteCarlo_Antithetic::Simulate(double start, double end, size_t steps) 
{

	std::cout << "American MC with antithetic" << std::endl;
	matrix weights = Payoff->GetWeights();
	//weights.Print();
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
	matrix paths_anti;


	for (size_t s = 0; s < m_Simulation; s++)
	{

		//std::cout << "simulation " << s << std::endl;
		if (s % 2 == 0)
		{
			//std::cout << "normal" << std::endl;
			x_diffusion->Simulate(start, end, steps);
			paths = x_diffusion->GetAllPaths();
			//std::cout << paths.nb_cols() << std::endl;
			//std::cout << paths.nb_rows() << std::endl;
		}
		else
		{
			//std::cout << "anti" << std::endl;
			paths_anti = x_diffusion->GetAllPathsAnti();
			//paths_anti.Print();
			//std::cout << paths_anti.nb_cols() << std::endl;
			//std::cout << paths_anti.nb_rows() << std::endl;
		}
		//paths.Print();
		chm = weights * paths;

		//chm.Print();
		//Index.addrows(chm.GetMatrix());
		//appendrow(Index, chm);//void function that modifies the Index matrix 
		//Index.Print();

		for (size_t i = 0; i < steps; i++)
		{
			Index(s, i) = chm(0, i);

		}

		simulated_price(s, 0) = Payoff->operator()(Index(s, Index.nb_cols() - 1));
		//simulated_price.Print();
		//populate the vector at maturity of the payoff

		//std::cout << "here" << std::endl;
	}

	// Regression part 

	//Index.Print();

	for (size_t t = steps - 2; t > 0; t--)
	{


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

		//2) construct Phi that contains the polynominial transformation at each order 

		//for (size_t i = 0; i < Index.nb_rows(); i++) 
		//{
		//	
		//	
		//			
		//}
		//It = Index.area(Index.nb_rows() - 1, t, { 0,t });
		//std::cout << "It" << std::endl;
		//It.Print();
		for (size_t b = 0; b < Phi.size(); b++)

		{
			poly_order_n = Phi[b]->operator()(It);
			//std::cout << "plynome " << b << std::endl;
			//poly_order_n.Print();

			for (size_t i = 0; i < poly_order_n.nb_rows(); i++)
			{
				Phit(i, b) = poly_order_n(i, 0);
			}


		}

		//std::cout << "Phit " << std::endl;
		//Phit.Print();

		Phit_T = transpose(Phit);

		SDP_Phi = Phit_T * Phit;

		Phi_V = Phit_T * V;

		inv_SDP_Phi = Inverse(SDP_Phi, SDP_Phi.nb_cols());

		Beta_hat = inv_SDP_Phi * Phi_V;

		C_hat = Phit * Beta_hat;

		for (size_t i = 0; i < C_hat.nb_rows(); i++)
		{

			if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(i, t)) > C_hat(i, 0)))
			{
				simulated_price(i, 0) = Payoff->operator()(Index(i, t));
				early_exec.push_back(i);
			}

			{


				simulated_price(i, 0) = df * simulated_price(i, 0);
			}


		}
	}

	MC_price = simulated_price.mean();

	MC_variance = simulated_price.variance();

};
>>>>>>> LSM working in vanilla, CV and antithetic
