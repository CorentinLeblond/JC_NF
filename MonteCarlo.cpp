#include "MonteCarlo.h"
#include <iostream>

//MC file

/*This file contains all algorithms used to do Monte Carlo simulations, with European or American payoff */

//Getters
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

//Optimal number of simulation with confidence level 99%
void MonteCarlo::OptimalNbSimul(const double& errortolerated)
{
	//The two-sided 99% quantile is:
	double quantile = 2.57;
	
	//It returns the minimum number of simulations required to make the estimated price below the error tolerated with confidence level 99%
	m_Simulation = (size_t) ((quantile*quantile*MC_variance)/(errortolerated*errortolerated));
	//The process is tho run a first companion MC on a few simulations, then estimate the variance, then define the minimul number and run the MC over this computed number
};

/////////////////// EU ///////////////////////////////////////////////////////

/* The constructor takes a number of simulation, a payoff, and a diffusion as inputs */
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
//European Basket Options
/////////////////////////////////////////////////////////////

EuropeanBasket::EuropeanBasket(size_t nbSimu, PayOffBasket* inputPayoff, RandomProcess* diffusion):
	MonteCarloEuropean(nbSimu,inputPayoff,diffusion)
{
};


void EuropeanBasket::Simulate(double start, double end, size_t steps)
{
	std::cout << "MC European Basket" << std::endl;
	simulated_price.Clear(); //Clear the simulated price matrix
	simulated_price.Resize(m_Simulation,1); //To store all prices simulated
	//Loop for all simulated trajectories
	for(size_t s = 0; s < m_Simulation; ++s) 
	{
		//Generate random process for N correlated assets 
		m_diffusion->Simulate(start, end, steps);
		//And store it
		matrix paths = m_diffusion->GetAllPaths();
		
		//Used to store underlying prices at maturity
		matrix maturity_spot(paths.nb_rows(), 1);
		
		//Loop to store underlying prices at maturity
		for (size_t i = 0; i < paths.nb_rows(); i++) 
		{
			//get the vector at maturity 
			maturity_spot(i, 0) = paths(i, paths.nb_cols() - 1);
			
		}
		//Compute the payoff
		simulated_price(s,0) = Payoff->operator()(maturity_spot);
	}

	//After all simulations being done, we store the average payoff discounted
	MC_price = exp(-r*end)*simulated_price.mean();
	//We also store the variance of the simulation
	MC_variance = simulated_price.variance();

};

//MC simulation with a pseudo-control variate 
EuropeanBasket_controlvariable::EuropeanBasket_controlvariable(size_t nbSimu, PayOffBasket* inputPayoff, PayOffBasket* Payoff_control,
																RandomProcess* diffusion,double inputclosedPrice)
	:
	MonteCarloEuropean(nbSimu,inputPayoff,diffusion),
	CPayoff(Payoff_control), //Payoff object for the control variate
	ExpPriceClsForm(inputclosedPrice) //Price of the control variate payoff computed with closed formula
{
};

void EuropeanBasket_controlvariable::Simulate(double start, double end, size_t steps)
{
	std::cout << "MC European Basket with the CV" << std::endl;
	simulated_price.Clear(); //Clear the simulated price matrix
	simulated_price.Resize(m_Simulation,1); //To store all prices simulated
	
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
		
		//Here is applied the CV mechanism, we substract the CV payoff, which is always below our basket call payoff, and add the closed formula price
		simulated_price(s, 0) = Payoff->operator()(maturity_spot)- CPayoff->operator()(maturity_spot)+ExpPriceClsForm;

	}

	MC_price = exp(-r*end)*simulated_price.mean();
	MC_variance = simulated_price.variance();
	
};

/* To make the antithetic class work well, the RandomProcess* input must a corresponding BS N dimensions Antithetic object
We use the fact that this object store each time a path and its antithetic in two different variables that are retreived 
thanks to two methods: getallpaths and getallpathsanti.
 */
EuropeanBasket_Antithetic::EuropeanBasket_Antithetic(size_t nbSimu, PayOffBasket* inputPayoff, RandomProcess* diffusion):
	MonteCarloEuropean(nbSimu,inputPayoff,diffusion)
{
};

void EuropeanBasket_Antithetic::Simulate(double start, double end, size_t steps)
{
	std::cout << "MC European Basket with Antithetic"<< std::endl;
	
	simulated_price.Clear(); //Clear the simulated price matrix
	simulated_price.Resize(m_Simulation/2,1);//To store all prices simulated
	
	simulated_price_Anti.Clear(); //Clear the simulated price matrix
	simulated_price_Anti.Resize(m_Simulation/2,1);
	
	//The average price matrix contains the average between a price for a given path, and the antithetic price for the corresponding antithetic path
	average_price.Clear(); //Clear the average price matrix,
	average_price.Resize(m_Simulation/2,1);
	
	size_t k = 0;

	for (size_t s = 0; s < m_Simulation; s++) 
	{
		//We start from s =0, so the first loop simulates a path, the second one only uses its antithetic path, and then again we simulate a path and its antithetic etc.
		if(s % 2 ==0)
		{
			m_diffusion->Simulate(start, end, steps);
			paths = m_diffusion->GetAllPaths();
		}
		else
		{
			paths = m_diffusion->GetAllPathsAnti();
		}

		matrix maturity_spot(paths.nb_rows(), 1);

		for (size_t i = 0; i < paths.nb_rows(); i++) 
		{			//std::cout << paths(i, paths.nb_col:endl;
			maturity_spot(i, 0) = paths(i, paths.nb_cols() - 1);
			//get the vector at maturity 
		}

		//Same mechanism, we store in two different matrices the prices and do an average once we have the price and the antithetic one.
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
														RandomProcess* diffusion,double inputclosedPrice)
	:
	EuropeanBasket_Antithetic(nbSimu,Payoff,diffusion),
	CPayoff(Payoff_control),
	ExpPriceClsForm(inputclosedPrice)
{
};

void EuropeanBasket_Antithetic_CV::Simulate(double start, double end, size_t steps)
{

	std::cout << "MC European Basket with Anti and CV"<< std::endl;
	
	simulated_price.Clear(); //Clear the simulated price matrix
	simulated_price.Resize(m_Simulation/2,1);//To store all prices simulated
	
	simulated_price_Anti.Clear(); //Clear the simulated price matrix
	simulated_price_Anti.Resize(m_Simulation/2,1);
	
	//The average price matrix contains the average between a price for a given path, and the antithetic price for the corresponding antithetic path
	average_price.Clear(); //Clear the average price matrix,
	average_price.Resize(m_Simulation/2,1);
	
	size_t k = 0;
	
	for (size_t s = 0; s < m_Simulation; s++) 
	{
		//std::cout << "simulation " << s << std::endl;
		if(s % 2 ==0)
		{
			m_diffusion->Simulate(start, end, steps);
			paths = m_diffusion->GetAllPaths();
		}
		else
		{
			paths = m_diffusion->GetAllPathsAnti();
		}
	
		matrix maturity_spot(paths.nb_rows(), 1);

		for (size_t i = 0; i < paths.nb_rows(); i++) 
		{
			maturity_spot(i, 0) = paths(i, paths.nb_cols() - 1);			
			//get the vector at maturity 
		}

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

	//Create the matrix at each time step of the polynome transformation in order to create the regressors
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


	Beta_hat = inv_SDP_Phi * Phi_V;

	C_hat = Phit * Beta_hat;

	Phit.Clear();
	poly_order_n.Clear();
	Phit_T.Clear();
	SDP_Phi.Clear();
	Phi_V.Clear();
	inv_SDP_Phi.Clear();
	Beta_hat.Clear();
	
	return C_hat;

};


matrix AmericanMonteCarlo::GetEarlyExec()
{
	//return the matrix with optimal stopping time at each time step 
	matrix res(stopping_time.size(), 1);

	for (size_t j = 0; j < stopping_time.size(); j++)
	{
		res(j, 0) = stopping_time[j];
	};

	return res;

};

matrix AmericanMonteCarlo::Final_Discounting(matrix simulated_price, std::vector<size_t> optimal_stopping_time, double dt, double rate) 
{
	//discount each simulated price at the appropriate discount factor taking into account the optimal stopping time 
	matrix final_discounted_payoff(simulated_price.nb_rows(), 1);

	for (size_t i = 0; i < simulated_price.nb_rows(); i++)
	{
		final_discounted_payoff(i, 0) = exp(-rate * optimal_stopping_time[i] * dt) * simulated_price(i, 0);
	}

	return final_discounted_payoff;

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
	matrix C_hat;
	matrix chm(1,steps);
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
		
		simulated_price(s, 0) = Payoff->operator()(Index(s,Index.nb_cols()-1)); //exp(-r*(end-start))*
		//populate the vector at maturity of the payoff

	}

	stopping_time.assign(m_Simulation, steps); // initialization of the stopping time at maturity for all paths

	double dt_sde = m_diffusion->Get_Dt();
	double df = exp(-r*dt_sde);

	// Regression part 
	for (size_t t = steps - 2; t > 0; t--)
	{
		std::vector<double> It;
		std::vector<double> V;
		std::vector<size_t> is_ITM;


		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			

			if (Payoff->operator()(Index(i, t)) > 0)
			{
				ITM(i, 0) = 1;
				It.push_back(Index(i, t)); //start as a vector where each time we add only the rows that are ITM
				V.push_back(df * simulated_price(i, 0));
				//V.push_back(df * Payoff->operator()(Index(i,t)));
				is_ITM.push_back(i);
			}
			else { ITM(i, 0) = -1; };
		}

		if (It.size() == 0) {
			std::cout << "Warning 0 paths ITM" << std::endl;
			goto increment;
		}
		{
		matrix It_mat(It, It.size()); //by construction this will be a column vector
		matrix V_mat(V, V.size());	//by construction this will be a column vector

		C_hat = C_Hat_regression(It_mat, V_mat);

		size_t c = 0;

		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{

			if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(is_ITM[c], t)) > C_hat(c, 0)))
			{

				simulated_price(i, 0) = Payoff->operator()(Index(is_ITM[c], t)); 
				stopping_time[i] = t;
				c += 1;
			}

		

			{
				simulated_price(i, 0) = simulated_price(i, 0);
				stopping_time[i] = stopping_time[i];
			}

		}

		It_mat.Clear();
		V_mat.Clear();
		It.clear();
		V.clear();
		is_ITM.clear();
	}
	increment:
		//increment is to take into account the possible fact that there are no path ITM and avoid creating empty matrix 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			simulated_price(i, 0) = simulated_price(i, 0);
		}

	}

	//After the backward regression simulated_price contains the cash_flow at the optimal stopping time 
	//stopping time contains the optimal stopping_time for each simulation
	//the function assigns to simulated_price the discounted payoff from the optimal stopping time to time 0
	simulated_price = Final_Discounting(simulated_price, stopping_time, dt_sde, r);
	

	MC_price = simulated_price.mean();

	MC_variance = simulated_price.variance();

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
	matrix C_hat;
	matrix chm(1, steps);
	matrix paths;


	matrix log_Index(m_Simulation, steps);
	matrix ITM2(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat_2;
	matrix chm_2(1, steps);
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
		simulated_price_CP(s,0) = CPayoff->operator()(log_Index(s, log_Index.nb_cols() - 1));
		//populate the vector at maturity of the payoff

	}

	stopping_time.assign(m_Simulation, steps); // initialization of the stopping time at maturity for all paths
	stopping_time_CP.assign(m_Simulation, steps);

	double dt_sde = m_diffusion->Get_Dt();
	double df = exp(-r * dt_sde);
	// Regression part 

	for (size_t t = steps - 2; t > 0; t--)
	{

		std::vector<double> V;
		std::vector<double> V_2;
		std::vector<size_t> is_ITM_2;
		std::vector<size_t> is_ITM;
		std::vector<double> It;
		std::vector<double> It_2;

		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			if (Payoff->operator()(Index(i, t)) > 0)
			{
				ITM(i, 0) = 1;
				It.push_back(Index(i, t));
				V.push_back(df * simulated_price(i, 0));
				is_ITM.push_back(i);
			}
			else { ITM(i, 0) = -1; };

			if (CPayoff->operator()(log_Index(i, t)) > 0)
			{
				ITM2(i, 0) = 1;
				It_2.push_back(log_Index(i, t));
				V_2.push_back(df * simulated_price_CP(i, 0));
				is_ITM_2.push_back(i);
			}
			else { ITM2(i, 0) = -1; };
		}

		if ((It.size() == 0)&&(It_2.size()==0))
		{
			std::cout << "Warning 0 paths ITM" << std::endl;
			goto  increment;
		}


		{
			matrix It_mat(It, It.size());
			matrix V_mat(V, V.size());
			matrix It_2_mat(It_2, It_2.size());
			matrix V_2_mat(V_2, V_2.size());
			C_hat = C_Hat_regression(It_mat, V_mat);
			C_hat_2 = C_Hat_regression(It_2_mat, V_2_mat);
			size_t c = 0;
			size_t c2 = 0;
			for (size_t i = 0; i < m_Simulation; i++)
			{

				if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(is_ITM[c], t)) > C_hat(c, 0)))
				{
					simulated_price(i, 0) = Payoff->operator()(Index(is_ITM[c], t));
					stopping_time[i] = t;
					c += 1;
				}

				{
					simulated_price(i, 0) = simulated_price(i, 0);
					stopping_time[i] = stopping_time[i];
				}

				if ((ITM2(i, 0) == 1) && (CPayoff->operator()(log_Index(is_ITM_2[c2], t)) > C_hat_2(c2, 0)))
				{
					simulated_price_CP(i, 0) = CPayoff->operator()(log_Index(is_ITM_2[c2], t));
					stopping_time_CP[i] = t;

					c2 += 1;
				}

				{

					simulated_price_CP(i, 0) = simulated_price_CP(i, 0);
					stopping_time_CP[i] = stopping_time_CP[i];
				}


			}

			

			It_mat.Clear();
			V_mat.Clear();
			It.clear();
			V.clear();
			is_ITM.clear();

			It_2_mat.Clear();
			V_2_mat.Clear();
			It_2.clear();
			V_2.clear();
			is_ITM_2.clear();
		}
	increment:
		for (size_t i = 0; i < m_Simulation; i++)
		{ 
			simulated_price(i, 0) = simulated_price(i, 0); 
			simulated_price_CP(i, 0) = simulated_price_CP(i, 0);
		}
	}

	simulated_price = Final_Discounting(simulated_price, stopping_time, dt_sde, r);
	simulated_price_CP = Final_Discounting(simulated_price_CP, stopping_time_CP, dt_sde, r);

	for (size_t i = 0; i < m_Simulation; i++)
	{
		CV(i, 0) = simulated_price(i, 0) - simulated_price_CP(i, 0) + ExpPriceClsForm;
	}

	

	MC_price = CV.mean();

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
	matrix C_hat;
	matrix chm(1, steps);
	matrix paths;
	matrix paths_anti;


	matrix chm_anti(1, steps);
	matrix Index_anti(m_Simulation/2, steps);
	matrix ITM_anti(m_Simulation/2, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat_anti;

	size_t k_a = 0;

	for (size_t s = 0; s < m_Simulation; s++)
	{

		if (s % 2 == 0)
		{
		m_diffusion->Simulate(start, end, steps);
			paths = m_diffusion->GetAllPaths();
			chm = weights * paths;

			for (size_t i = 0; i < steps; i++)
			{
				Index(k_a, i) = chm(0, i);

			}

			simulated_price(k_a, 0) = Payoff->operator()(Index(k_a, Index.nb_cols() - 1));


			

		}
		else
		{
			paths_anti = m_diffusion->GetAllPathsAnti();
			chm_anti = weights * paths_anti;


			for (size_t i = 0; i < steps; i++)
			{
				Index_anti(k_a, i) = chm_anti(0, i);

			}

			simulated_price_anti(k_a, 0) = Payoff->operator()(Index_anti(k_a, Index_anti.nb_cols() - 1));
			k_a += 1;


		}
	}

	stopping_time.assign(m_Simulation/2, steps); // initialization of the stopping time at maturity for all paths
	stopping_time_anti.assign(m_Simulation/2, steps);
	double dt_sde = m_diffusion->Get_Dt();
	double df = exp(-r * dt_sde);
	// Regression part 

	for (size_t t = steps - 2; t > 0; t--)
	{

		std::vector<double> V;
		std::vector<double> V_anti;
		std::vector<size_t> is_ITM_anti;
		std::vector<size_t> is_ITM;
		std::vector<double> It;
		std::vector<double> It_anti;

		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
	

			if (Payoff->operator()(Index(i, t)) > 0)
			{
				ITM(i, 0) = 1;
				V.push_back(df * simulated_price(i, 0));
				It.push_back(Index(i, t));
				is_ITM.push_back(i);

			}
			else { ITM(i, 0) = -1; };

			if (Payoff->operator()(Index_anti(i, t)) > 0)
			{
				ITM_anti(i, 0) = 1;
				V_anti.push_back(df * simulated_price_anti(i, 0));
				It_anti.push_back(Index_anti(i, t));
				is_ITM_anti.push_back(i);
			}
			else { ITM_anti(i, 0) = -1; };
		}

		if ((It.size() == 0)&& (It_anti.size() == 0))
		{
			std::cout << "Warning 0 paths ITM" << std::endl;
			goto  increment;
		}

		{
			matrix It_mat(It, It.size());
			matrix It_mat_anti(It_anti, It_anti.size());

			matrix V_mat(V, V.size());
			matrix V_mat_anti(V_anti, V_anti.size());

			C_hat = C_Hat_regression(It_mat, V_mat);
			C_hat_anti = C_Hat_regression(It_mat_anti, V_mat_anti);

			size_t c = 0;
			size_t ca = 0;

			for (size_t i = 0; i < ITM.nb_rows(); i++)
			{
				if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(is_ITM[c], t)) > C_hat(c, 0)))
				{
					simulated_price(i, 0) = Payoff->operator()(Index(is_ITM[c], t));
					stopping_time[i] = t;
					c += 1;
				}

				{
					simulated_price(i, 0) = simulated_price(i, 0);
					stopping_time[i] = stopping_time[i];
				}

		

				if ((ITM_anti(i, 0) == 1) && (Payoff->operator()(Index_anti(is_ITM_anti[ca], t)) > C_hat_anti(ca, 0)))
				{
					simulated_price_anti(i, 0) = Payoff->operator()(Index_anti(is_ITM_anti[ca], t));
					stopping_time_anti[i] = t;
					ca += 1;
				}

				{
					simulated_price_anti(i, 0) = simulated_price_anti(i, 0);
					stopping_time_anti[i] = stopping_time_anti[i];
				}

			}

			It_mat.Clear();
			V_mat.Clear();
			It.clear();
			V.clear();
			is_ITM.clear();

			It_mat_anti.Clear();
			V_mat_anti.Clear();
			It_anti.clear();
			V_anti.clear();
			is_ITM_anti.clear();
		}
	increment:
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			simulated_price(i, 0) = simulated_price(i, 0);
			simulated_price_anti(i, 0) = simulated_price_anti(i, 0);
		}
	}

	simulated_price = Final_Discounting(simulated_price, stopping_time, dt_sde, r);
	simulated_price_anti = Final_Discounting(simulated_price_anti, stopping_time_anti, dt_sde, r);

	for (size_t i = 0; i < simulated_price.nb_rows(); i++) 
	{
		average_price(i, 0) = (simulated_price(i, 0) + simulated_price_anti(i, 0)) * 0.5;
		
	};

	MC_price = average_price.mean();

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
	matrix C_hat;
	matrix chm(1, steps);
	matrix paths;
	matrix paths_anti;


	matrix chm_anti(1, steps);
	matrix Index_anti(m_Simulation / 2, steps);
	matrix ITM_anti(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat_anti;



	matrix log_Index(m_Simulation/2, steps);
	matrix ITM2(m_Simulation/2, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat_2;
	matrix chm_2(1, steps);
	matrix log_spot;


	matrix log_Index_anti(m_Simulation/2, steps);
	matrix ITM2_anti(m_Simulation/2, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat_2_anti;
	matrix chm_2_anti(1, steps);
	matrix log_spot_anti;
	size_t k = 0;
	size_t k_a = 0;

	for (size_t s = 0; s < m_Simulation; s++)
	{

		//std::cout << "simulation " << s << std::endl;
		if (s % 2 == 0)
		{
			m_diffusion->Simulate(start, end, steps);
			paths = m_diffusion->GetAllPaths();
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
				Index(k_a, i) = chm(0, i);
				log_Index(k_a, i) = chm_2(0, i);

			}

			simulated_price(k_a, 0) = Payoff->operator()(Index(k_a, Index.nb_cols() - 1));
			simulated_price_CP(k_a, 0) = CPayoff->operator()(log_Index(k_a, log_Index.nb_cols() - 1));


		}
		else
		{
			paths_anti = m_diffusion->GetAllPathsAnti();
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

	stopping_time.assign(m_Simulation/2, steps);
	stopping_time_CP.assign(m_Simulation / 2, steps);
	stopping_time_anti.assign(m_Simulation / 2, steps);
	stopping_time_anti_CP.assign(m_Simulation / 2, steps);

	double dt_sde = m_diffusion->Get_Dt();
	double df = exp(-r *dt_sde);

	// Regression part 

	for (size_t t = steps - 2; t > 0; t--)
	{

		std::vector<double> V;
		std::vector<double> V_anti;
		std::vector<size_t> is_ITM_anti;
		std::vector<size_t> is_ITM;
		std::vector<double> It;
		std::vector<double> It_anti;


		std::vector<double> V_2;
		std::vector<double> V_anti_2;
		std::vector<size_t> is_ITM_anti_2;
		std::vector<size_t> is_ITM_2;
		std::vector<double> It_2;
		std::vector<double> It_anti_2;

		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			

			if (Payoff->operator()(Index(i,t)) > 0)
			{
				ITM(i, 0) = 1;
				V.push_back(df* simulated_price(i, 0));
				It.push_back(Index(i, t));
				is_ITM.push_back(i);
			}
			else { ITM(i, 0) = -1; };

			if (CPayoff->operator()(log_Index(i,t)) > 0)
			{
				ITM2(i, 0) = 1;
				V_2.push_back(df* simulated_price_CP(i, 0));
				It_2.push_back(log_Index(i, t));
				is_ITM_2.push_back(i);

			}
			else { ITM2(i, 0) = -1; };


			if (Payoff->operator()(Index_anti(i,t)) > 0)
			{
				ITM_anti(i, 0) = 1;
				V_anti.push_back(df* simulated_price_anti(i, 0));
				It_anti.push_back(Index_anti(i, t));
				is_ITM_anti.push_back(i);
			}
			else { ITM_anti(i, 0) = -1; };


			if (CPayoff->operator()(log_Index_anti(i, t)) > 0)
			{
				ITM2_anti(i, 0) = 1;
				V_anti_2.push_back(df* simulated_price_CP_anti(i, 0));
				It_anti_2.push_back(log_Index_anti(i, t));
				is_ITM_anti_2.push_back(i);
			}
			else { ITM2_anti(i, 0) = -1; };

		}

		if ((It.size() == 0) && (It_anti.size() == 0)&& (It_2.size() == 0) && (It_anti_2.size() == 0))
		{
			std::cout << "Warning 0 paths ITM" << std::endl;
			goto  increment;
		}

		{
			matrix It_mat(It, It.size());
			matrix It_mat_anti(It_anti, It_anti.size());

			matrix It2_mat(It_2, It_2.size());
			matrix It2_mat_anti(It_anti_2, It_anti_2.size());

			matrix V2_mat(V_2, V_2.size());
			matrix V2_mat_anti(V_anti_2, V_anti_2.size());

			matrix V_mat(V, V.size());
			matrix V_mat_anti(V_anti, V_anti.size());

			C_hat = C_Hat_regression(It_mat, V_mat);
			C_hat_2 = C_Hat_regression(It2_mat, V2_mat);
			C_hat_anti = C_Hat_regression(It_mat_anti, V_mat_anti);
			C_hat_2_anti = C_Hat_regression(It2_mat_anti, V2_mat_anti);


			size_t c = 0;
			size_t c2 = 0;
			size_t ca = 0;
			size_t ca_2 = 0;

			for (size_t i = 0; i < ITM.nb_rows(); i++)
			{

				if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(is_ITM[c], t)) > C_hat(c, 0)))
				{
					simulated_price(i, 0) = Payoff->operator()(Index(is_ITM[c], t));
					stopping_time[i] = t;
					c += 1;
				}
				else
				{
					simulated_price(i, 0) = simulated_price(i, 0);
					stopping_time[i] = stopping_time[i];
				}

				if ((ITM2(i, 0) == 1) && (CPayoff->operator()(log_Index(is_ITM_2[c2], t)) > C_hat_2(c2, 0)))
				{
					simulated_price_CP(i, 0) =CPayoff->operator()(log_Index(is_ITM_2[c2], t));
					stopping_time_CP[i] = t;
					c2 += 1;
				}
				else
				{
					simulated_price_CP(i, 0) = simulated_price_CP(i, 0);
					stopping_time_CP[i] = stopping_time_CP[i];
				}


				if ((ITM_anti(i, 0) == 1) && (Payoff->operator()(Index_anti(is_ITM_anti[ca], t)) > C_hat_anti(ca, 0)))
				{
					simulated_price_anti(i, 0) =Payoff->operator()(Index_anti(is_ITM_anti[ca], t));
					stopping_time_anti[i] = t;
					ca += 1;
				}
				else
				{
					simulated_price_anti(i, 0) = simulated_price_anti(i, 0);
					stopping_time_anti[i] = stopping_time_anti[i];
				}


				if ((ITM2_anti(i, 0) == 1) && (CPayoff->operator()(log_Index_anti(is_ITM_anti_2[ca_2], t)) > C_hat_2_anti(ca_2, 0)))
				{
					simulated_price_CP_anti(i, 0) = CPayoff->operator()(log_Index_anti(is_ITM_anti_2[ca_2], t));
					stopping_time_anti_CP[i] = t;
					ca_2 += 1;

				}
				else
				{
					simulated_price_CP_anti(i, 0) = simulated_price_CP_anti(i, 0);
					stopping_time_anti_CP[i] = stopping_time_anti_CP[i];
				}

			}

			It_mat.Clear();
			V_mat.Clear();
			It.clear();
			V.clear();
			is_ITM.clear();

			It2_mat.Clear();
			V2_mat.Clear();
			It_anti.clear();
			V_anti.clear();
			is_ITM_anti.clear();

			It_mat_anti.Clear();
			V_mat_anti.Clear();
			It_anti.clear();
			V_anti.clear();
			is_ITM_anti.clear();

			It2_mat_anti.Clear();
			V2_mat_anti.Clear();
			It_anti_2.clear();
			V_anti_2.clear();
			is_ITM_anti_2.clear(); 
		}
		increment:
			for (size_t i = 0; i < ITM.nb_rows(); i++) 
			{
				simulated_price(i, 0) = simulated_price(i,0);
				simulated_price_CP(i, 0) = simulated_price_CP(i, 0);
				simulated_price_CP_anti(i, 0) =simulated_price_CP_anti(i, 0);
				simulated_price_anti(i, 0) =simulated_price_anti(i, 0);
			
			}
	}

	simulated_price = Final_Discounting(simulated_price, stopping_time, dt_sde, r);
	simulated_price_CP = Final_Discounting(simulated_price_CP, stopping_time_CP, dt_sde, r);
	simulated_price_anti = Final_Discounting(simulated_price_anti, stopping_time_anti, dt_sde, r);
	simulated_price_CP_anti = Final_Discounting(simulated_price_CP_anti, stopping_time_anti_CP, dt_sde, r);

	for (size_t i = 0; i < simulated_price.nb_rows(); i++)
	{
		temp = simulated_price(i, 0) - simulated_price_CP(i, 0);
		temp_anti = simulated_price_anti(i, 0) - simulated_price_CP_anti(i, 0);

		average_price(i, 0) = (temp + temp_anti) * 0.5 + ExpPriceClsForm;

		temp = 0.;
		temp_anti = 0;
	}



	MC_price = average_price.mean();

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


	std::cout << "Bermudean BasketOption LS" << std::endl;
	simulated_price.Resize(m_Simulation, 1);
	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation, exec_schedule.nb_rows());
	matrix ITM(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat;
	matrix chm(1, exec_schedule.nb_rows());
	matrix paths;
	matrix interpolated_paths;
	double dt_sde;
	double df;

	for (size_t s = 0; s < m_Simulation; s++)
	{
		//std::cout << "simulation " << s << std::endl;
		m_diffusion->Simulate(start, end, steps);
		paths = m_diffusion->GetAllPaths();
		interpolated_paths = wkday->index_executed(paths, exec_schedule, steps);

		
		//interpolated_paths.Print();
		chm = weights * interpolated_paths;

		for (size_t i = 0; i < chm.nb_cols(); i++)
		{
			Index(s, i) = chm(0, i);

		}
		simulated_price(s, 0) = Payoff->operator()(Index(s, Index.nb_cols() - 1));
		//populate the vector at maturity of the payoff

	}

	
	dt_sde = m_diffusion->Get_Dt();

	for (size_t j = 1; j < exec_schedule.nb_rows(); j++)
	{

		Dt_schedule(j - 1, 0) = (ceil(exec_schedule(j, 0) * steps) - ceil(exec_schedule(j - 1, 0) * steps)) * dt_sde;

	}


	size_t Isteps = Index.nb_cols() -2; //Index will have as much columns as there are execution dates.

	stopping_time.assign(m_Simulation, ceil(exec_schedule(exec_schedule.nb_rows()-1, 0) * steps));
	//the first available stoppping time is the last execution date available

	// Regression part
	for (size_t t = Isteps; t > 0; t--)
	{
		std::vector<double> It;
		std::vector<double> V;
		std::vector<size_t> is_ITM;

		df = exp(-r * Dt_schedule(t, 0));
		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			if (Payoff->operator()(Index(i, t)) > 0)
			{
				ITM(i, 0) = 1;
				It.push_back(Index(i,t));
				V.push_back(df * simulated_price(i, 0));
				is_ITM.push_back(i);
				
			}
			else { ITM(i, 0) = -1; };
		}
		if (It.size() == 0) 
		{
			std::cout << "Warning, 0 paths ITM" << std::endl; goto increment;
		};
		{
			matrix It_mat(It, It.size());
			matrix V_mat(V, V.size());

			C_hat = C_Hat_regression(It_mat, V_mat); //the function will apply the polynomes, create Beta_Hat, and return C_hat as a matrix

			size_t c = 0;
			for (size_t i = 0; i < ITM.nb_rows(); i++)
			{

				if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(is_ITM[c], t)) > C_hat(c, 0)))
				{
					simulated_price(i, 0) = Payoff->operator()(Index(is_ITM[c], t));
					stopping_time[i] = t;
					c += 1;
				}

				{
					simulated_price(i, 0) = simulated_price(i, 0);
					stopping_time[i] = stopping_time[i];
				}


			}

			It.clear();
			V.clear();
			It_mat.Clear();
			V_mat.Clear();
			is_ITM.clear();
		}
	increment: 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
			simulated_price(i, 0) = simulated_price(i, 0);
		}
	}

	simulated_price = Final_Discounting(simulated_price, stopping_time, dt_sde, r);

	MC_price = simulated_price.mean();

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
	

	

	std::cout << "Bermudean MC LS with Control Variable" << std::endl;

	simulated_price.Resize(m_Simulation, 1);
	simulated_price_CP.Resize(m_Simulation, 1);
	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation, exec_schedule.nb_rows());
	matrix ITM(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat;
	matrix chm;
	matrix paths;

	matrix log_Index(m_Simulation, exec_schedule.nb_rows());
	matrix ITM2(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat_2;
	matrix chm_2;
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
		interpolated_paths = wkday->index_executed(paths, exec_schedule, steps);

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
		simulated_price_CP(s, 0) = CPayoff->operator()(log_Index(s, log_Index.nb_cols() - 1));
		//populate the vector at maturity of the payoff

	}

	dt_sde = m_diffusion->Get_Dt();

	for (size_t j = 1; j < exec_schedule.nb_rows(); j++)
	{

		Dt_schedule(j - 1, 0) = (ceil(exec_schedule(j, 0) * steps) - ceil(exec_schedule(j - 1, 0) * steps)) * dt_sde;

	}

	double last_dt = ceil(exec_schedule(0, 0) * steps) * dt_sde;

	stopping_time.assign(m_Simulation, ceil(exec_schedule(exec_schedule.nb_rows() - 1, 0) * steps));
	stopping_time_CP.assign(m_Simulation, ceil(exec_schedule(exec_schedule.nb_rows() - 1, 0) * steps));

	// Regression part

	size_t Isteps = Index.nb_cols() - 2;

	for (size_t t = Isteps; t > 0; t--)
	{
		std::vector<double> V;
		std::vector<double> V_2;
		std::vector<size_t> is_ITM_2;
		std::vector<size_t> is_ITM;
		std::vector<double> It;
		std::vector<double> It_2;

		df = exp(-r * Dt_schedule(t, 0));
		// 1) separate ITM path from others 

		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{

			if (Payoff->operator()(Index(i, t)) > 0)
			{
				ITM(i, 0) = 1;
				V.push_back(df * simulated_price(i, 0));
				It.push_back(Index(i, t));
				is_ITM.push_back(i);
			}
			else { ITM(i, 0) = -1; };

			if (CPayoff->operator()(log_Index(i, t)) > 0)
			{
				ITM2(i, 0) = 1;
				V_2.push_back(df * simulated_price_CP(i, 0));
				It_2.push_back(log_Index(i, t));
				is_ITM_2.push_back(i);
			}
			else { ITM2(i, 0) = -1; };
		}

			if ((It.size() == 0) && (It_2.size() == 0))
			{
				std::cout << "Warning 0 paths ITM" << std::endl;
				goto increment;
			}

			{
				matrix It_mat(It, It.size());
				matrix V_mat(V, V.size());
				matrix It_mat_2(It_2, It_2.size());
				matrix V_mat_2(V_2, V_2.size());
				C_hat = C_Hat_regression(It_mat, V_mat);
				C_hat_2 = C_Hat_regression(It_mat_2, V_mat_2);

				size_t c = 0;
				size_t c2 = 0;

				for (size_t i = 0; i < ITM.nb_rows(); i++)
				{

					if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(is_ITM[c], t)) > C_hat(c, 0)))
					{
						simulated_price(i, 0) = Payoff->operator()(Index(is_ITM[c], t));
						stopping_time[i] = t;
						c += 1;
					}

					{
						simulated_price(i, 0) = simulated_price(i, 0);
						stopping_time[i] = stopping_time[i];
					}

			



					if ((ITM2(i, 0) == 1) && (CPayoff->operator()(log_Index(is_ITM_2[c2], t)) > C_hat_2(c2, 0)))
					{
						simulated_price_CP(i, 0) = CPayoff->operator()(log_Index(is_ITM_2[c2], t));
						stopping_time_CP[i] = t;
						c2 += 1;
					}

					{
						simulated_price_CP(i, 0) =simulated_price_CP(i, 0);
						stopping_time_CP[i] = stopping_time_CP[i];
					}

					


				}

				V_2.clear();
				V.clear();
				It.clear();
				It_2.clear();

				It_mat.Clear();
				V_mat.Clear();
				is_ITM.clear();
				is_ITM_2.clear();

				It_mat_2.Clear();
				V_mat_2.Clear();

			}

		increment:

			for (size_t i = 0; i < ITM.nb_rows(); i++)
			{

				simulated_price(i, 0) = simulated_price(i, 0);
				simulated_price_CP(i, 0) = simulated_price_CP(i, 0);
			}
		}
	
		simulated_price = Final_Discounting(simulated_price, stopping_time, dt_sde, r);
		simulated_price_CP = Final_Discounting(simulated_price_CP, stopping_time_CP, dt_sde, r);

		for (size_t i = 0; i < simulated_price.nb_rows(); i++)
		{

			CV(i, 0) = simulated_price(i, 0) - simulated_price_CP(i, 0) + ExpPriceClsForm;
		}


	MC_price = CV.mean();
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


	std::cout << "Bermudean MC with antithetic" << std::endl;

	simulated_price.Resize(m_Simulation / 2, 1);
	simulated_price_anti.Resize(m_Simulation / 2, 1);
	average_price.Resize(m_Simulation / 2, 1);

	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation / 2, exec_schedule.nb_rows());
	matrix ITM(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat;
	matrix chm;
	matrix paths;
	matrix interpolated_paths;
	matrix paths_anti;
	matrix interpolated_paths_anti;
	

	matrix chm_anti;
	matrix Index_anti(m_Simulation / 2, exec_schedule.nb_rows());
	matrix ITM_anti(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat_anti;
	double dt_sde;
	double df;
	size_t k_a = 0;

	for (size_t s = 0; s < m_Simulation; s++)
	{

		if (s % 2 == 0)
		{
			m_diffusion->Simulate(start, end, steps);
			paths = m_diffusion->GetAllPaths();
			interpolated_paths = wkday->index_executed(paths, exec_schedule, steps);
			chm = weights * interpolated_paths;

			for (size_t i = 0; i < chm.nb_cols(); i++)
			{
				Index(k_a, i) = chm(0, i);

			}

			simulated_price(k_a, 0) = Payoff->operator()(Index(k_a, Index.nb_cols() - 1));

		}
		else
		{
			paths_anti = m_diffusion->GetAllPathsAnti();
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

	dt_sde = m_diffusion->Get_Dt();

	for (size_t j = 1; j < exec_schedule.nb_rows(); j++)
	{

		Dt_schedule(j - 1, 0) = (ceil(exec_schedule(j, 0) * steps) - ceil(exec_schedule(j - 1, 0) * steps)) * dt_sde;

	}

	double last_dt = ceil(exec_schedule(0, 0) * steps) * dt_sde;

	stopping_time.assign(m_Simulation/2, ceil(exec_schedule(exec_schedule.nb_rows() - 1, 0) * steps));
	stopping_time_anti.assign(m_Simulation / 2, ceil(exec_schedule(exec_schedule.nb_rows() - 1, 0) * steps));

	// Regression part 
	size_t Isteps = Index.nb_cols() - 2;

	for (size_t t = Isteps; t > 0; t--)
	{

		std::vector<double> V;
		std::vector<double> V_anti;
		std::vector<size_t> is_ITM_anti;
		std::vector<size_t> is_ITM;
		std::vector<double> It;
		std::vector<double> It_anti;


		df = exp(-r * Dt_schedule(t, 0));

		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{

			if (Payoff->operator()(Index(i, t)) > 0)
			{
				ITM(i, 0) = 1;
				It.push_back(Index(i, t));
				V.push_back(df * simulated_price(i, 0));
				is_ITM.push_back(i);
			}
			else { ITM(i, 0) = -1; };

			if (Payoff->operator()(Index_anti(i, t)) > 0)
			{
				ITM_anti(i, 0) = 1;
				It_anti.push_back(Index_anti(i, t));
				V_anti.push_back(df * simulated_price_anti(i, 0));
				is_ITM_anti.push_back(i);
			}
			else { ITM_anti(i, 0) = -1; };
		}

		if ((It.size() == 0) && (It_anti.size() == 0))
		{
			std::cout << "Warning 0 paths ITM" << std::endl;
			goto increment;
		}

		{
			matrix It_mat(It, It.size());
			matrix V_mat(V, V.size());
			matrix It_mat_anti(It_anti, It_anti.size());
			matrix V_mat_anti(V_anti, V_anti.size());

			C_hat = C_Hat_regression(It_mat, V_mat);

			C_hat_anti = C_Hat_regression(It_mat_anti, V_mat_anti);

			size_t c = 0;
			size_t ca = 0;

			for (size_t i = 0; i < ITM.nb_rows(); i++)
			{


				if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(is_ITM[c], t)) > C_hat(c, 0)))
				{
					simulated_price(i, 0) = Payoff->operator()(Index(is_ITM[c], t));
					stopping_time[i] = t;
					c += 1;
				}

				{
					simulated_price(i, 0) = simulated_price(i, 0);
					stopping_time[i] = stopping_time[i];
				}


				if ((ITM_anti(i, 0) == 1) && (Payoff->operator()(Index_anti(is_ITM_anti[ca], t)) > C_hat_anti(ca, 0)))
				{
					simulated_price_anti(i, 0) = Payoff->operator()(Index_anti(is_ITM_anti[ca], t));
					stopping_time_anti[i] = t;
					ca += 1;
				}

				{

					simulated_price_anti(i, 0) = simulated_price_anti(i, 0);
					stopping_time_anti[i] = stopping_time_anti[i];
				}

				
			}

			V_anti.clear();
			V.clear();
			It.clear();
			It_anti.clear();

			It_mat.Clear();
			V_mat.Clear();
			is_ITM.clear();
			is_ITM_anti.clear();

			It_mat_anti.Clear();
			V_mat_anti.Clear();
		}

	increment: 

		for (size_t i = 0; i < ITM.nb_rows(); i++) 
		{
			simulated_price(i, 0) = simulated_price(i, 0);
			simulated_price_anti(i, 0) = simulated_price_anti(i, 0);
			
		}
		
	}

	simulated_price = Final_Discounting(simulated_price, stopping_time, dt_sde, r);
	simulated_price_anti = Final_Discounting(simulated_price_anti, stopping_time_anti, dt_sde, r);

	for (size_t i = 0; i < simulated_price.nb_rows(); i++)
	{

		average_price(i, 0) = (simulated_price(i, 0) + simulated_price_anti(i, 0)) * 0.5;
	}


	
	MC_price = average_price.mean();

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


	matrix weights = Payoff->GetWeights();
	matrix Index(m_Simulation / 2, exec_schedule.nb_rows());
	matrix ITM(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat;
	matrix chm;
	matrix paths;
	matrix paths_anti;


	matrix chm_anti;
	matrix Index_anti(m_Simulation / 2, exec_schedule.nb_rows());
	matrix ITM_anti(m_Simulation / 2, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat_anti;



	matrix log_Index(m_Simulation, exec_schedule.nb_rows());
	matrix ITM2(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat_2;
	matrix chm_2;
	matrix log_spot;


	matrix log_Index_anti(m_Simulation, exec_schedule.nb_rows());
	matrix ITM2_anti(m_Simulation, 1); //need to separate which paths are ITM at each iteration 
	matrix C_hat_2_anti;
	matrix chm_2_anti;
	matrix log_spot_anti;

	matrix interpolated_paths;
	matrix interpolated_paths_anti;
	double dt_sde;
	double df;

	size_t k = 0;
	size_t k_a = 0;

	for (size_t s = 0; s < m_Simulation; s++)
	{

		
		if (s % 2 == 0)
		{
			m_diffusion->Simulate(start, end, steps);
			paths = m_diffusion->GetAllPaths();
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

	dt_sde = m_diffusion->Get_Dt();

	for (size_t j = 1; j < exec_schedule.nb_rows(); j++)
	{

		Dt_schedule(j - 1, 0) = (ceil(exec_schedule(j, 0)*steps) - ceil(exec_schedule(j - 1, 0)*steps)) * dt_sde;

	}

	double last_dt = ceil(exec_schedule(0, 0)) * dt_sde;
	stopping_time.assign(m_Simulation/2, ceil(exec_schedule(exec_schedule.nb_rows() - 1, 0) * steps));
	stopping_time_CP.assign(m_Simulation / 2, ceil(exec_schedule(exec_schedule.nb_rows() - 1, 0) * steps));
	stopping_time_anti.assign(m_Simulation / 2, ceil(exec_schedule(exec_schedule.nb_rows() - 1, 0) * steps));
	stopping_time_CP_anti.assign(m_Simulation / 2, ceil(exec_schedule(exec_schedule.nb_rows() - 1, 0) * steps));
	// Regression part 

	size_t Isteps = Index.nb_cols() - 2;

	for (size_t t = Isteps; t > 0; t--)
	{

		std::vector<double> V;
		std::vector<double> V_anti;
		std::vector<size_t> is_ITM_anti;
		std::vector<size_t> is_ITM;
		std::vector<double> It;
		std::vector<double> It_anti;


		std::vector<double> V_2;
		std::vector<double> V_anti_2;
		std::vector<size_t> is_ITM_anti_2;
		std::vector<size_t> is_ITM_2;
		std::vector<double> It_2;
		std::vector<double> It_anti_2;


		df = exp(-r * Dt_schedule(t, 0));
		// 1) separate ITM path from others 
		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
	

			if (Payoff->operator()(Index(i, t)) > 0)
			{
				ITM(i, 0) = 1;
				V.push_back(df * simulated_price(i, 0));
				It.push_back(Index(i, t));
				is_ITM.push_back(i);


			}
			else { ITM(i, 0) = -1; };

			if (CPayoff->operator()(log_Index(i, t)) > 0)
			{
				ITM2(i, 0) = 1;
				V_2.push_back(df * simulated_price_CP(i, 0));
				It_2.push_back(log_Index(i, t));
				is_ITM_2.push_back(i);

			}
			else { ITM2(i, 0) = -1; };


			if (Payoff->operator()(Index_anti(i, t)) > 0)
			{
				ITM_anti(i, 0) = 1;
				//V_anti.push_back(df * simulated_price_anti(i, 0));
				V_anti.push_back(df* Payoff->operator()(Index_anti(i, t)));
				It_anti.push_back(Index_anti(i, t));
				is_ITM_anti.push_back(i);

			}
			else { ITM_anti(i, 0) = -1; };


			if (CPayoff->operator()(log_Index_anti(i, t)) > 0)
			{
				ITM2_anti(i, 0) = 1;
				V_anti_2.push_back(df * simulated_price_CP_anti(i, 0));
				It_anti_2.push_back(log_Index_anti(i, t));
				is_ITM_anti_2.push_back(i);

			}
			else { ITM2_anti(i, 0) = -1; };

		}

		if ((It.size() == 0) && (It_2.size() == 0) && (It_anti.size() == 0) && (It_anti_2.size() == 0))
		{
			std::cout << "Warning 0 paths ITM" << std::endl;
			goto increment;
		}


		{
			matrix It_mat(It, It.size());
			matrix It_mat_anti(It_anti, It_anti.size());

			matrix It2_mat(It_2, It_2.size());
			matrix It2_mat_anti(It_anti_2, It_anti_2.size());

			matrix V2_mat(V_2, V_2.size());
			matrix V2_mat_anti(V_anti_2, V_anti_2.size());

			matrix V_mat(V, V.size());
			matrix V_mat_anti(V_anti, V_anti.size());

			C_hat = C_Hat_regression(It_mat, V_mat);
			C_hat_2 = C_Hat_regression(It2_mat, V2_mat);
			C_hat_anti = C_Hat_regression(It_mat_anti, V_mat_anti);
			C_hat_2_anti = C_Hat_regression(It2_mat_anti, V2_mat_anti);

			size_t c = 0;
			size_t c2 = 0;
			size_t ca = 0;
			size_t ca_2 = 0;

			for (size_t i = 0; i < ITM.nb_rows(); i++)
			{

				if ((ITM(i, 0) == 1) && (Payoff->operator()(Index(is_ITM[c], t)) > C_hat(c, 0)))
				{
					simulated_price(i, 0) = Payoff->operator()(Index(is_ITM[c], t));
					stopping_time[i] = t;
					c += 1;
				}
				else
				{

					simulated_price(i, 0) = simulated_price(i, 0);
					stopping_time[i] = stopping_time[i];
				}

				if ((ITM2(i, 0) == 1) && (CPayoff->operator()(log_Index(is_ITM_2[c2], t)) > C_hat_2(c2, 0)))
				{
					simulated_price_CP(i, 0) = CPayoff->operator()(log_Index(is_ITM_2[c2], t));
					stopping_time_CP[i] = t;
					c2 += 1;
				}
				else
				{

					simulated_price_CP(i, 0) = simulated_price_CP(i, 0);
					stopping_time_CP[i] = stopping_time_CP[i];
				}


				if ((ITM_anti(i, 0) == 1) && (Payoff->operator()(Index_anti(is_ITM_anti[ca], t)) > C_hat_anti(ca, 0)))
				{
					simulated_price_anti(i, 0) = Payoff->operator()(Index_anti(is_ITM_anti[ca], t));
					stopping_time_anti[i] = t;
					ca += 1;
				}
				else
				{

					simulated_price_anti(i, 0) = simulated_price_anti(i, 0);
					stopping_time_anti[i] = stopping_time_anti[i];
				}


				if ((ITM2_anti(i, 0) == 1) && (CPayoff->operator()(log_Index_anti(is_ITM_anti_2[ca_2], t)) > C_hat_2_anti(ca_2, 0)))
				{
					simulated_price_CP_anti(i, 0) = CPayoff->operator()(log_Index_anti(is_ITM_anti_2[ca_2], t));
					stopping_time_CP_anti[i] = t;
					ca_2 += 1;
				}
				else
				{

					simulated_price_CP_anti(i, 0) = simulated_price_CP_anti(i, 0);
					stopping_time_CP_anti[i] = stopping_time_CP_anti[i];
				}


			}
		
			It_mat.Clear();
			V_mat.Clear();
			It.clear();
			V.clear();
			is_ITM.clear();

			It2_mat.Clear();
			V2_mat.Clear();
			It_anti.clear();
			V_anti.clear();
			is_ITM_anti.clear();

			It_mat_anti.Clear();
			V_mat_anti.Clear();
			It_anti.clear();
			V_anti.clear();
			is_ITM_anti.clear();

			It2_mat_anti.Clear();
			V2_mat_anti.Clear();
			It_anti_2.clear();
			V_anti_2.clear();
			is_ITM_anti_2.clear();
		
		}

	increment:

		for (size_t i = 0; i < ITM.nb_rows(); i++)
		{
				simulated_price(i, 0) = simulated_price(i,0);
				simulated_price_CP(i, 0) = simulated_price_CP(i, 0);

				simulated_price_anti(i, 0) = simulated_price_anti(i, 0);
				simulated_price_CP_anti(i, 0) = simulated_price_CP_anti(i, 0);
				
			}


	}

	simulated_price = Final_Discounting(simulated_price, stopping_time, dt_sde, r);
	simulated_price_anti = Final_Discounting(simulated_price_anti, stopping_time_anti, dt_sde, r);
	simulated_price_CP = Final_Discounting(simulated_price_CP, stopping_time_CP, dt_sde, r);
	simulated_price_CP_anti = Final_Discounting(simulated_price_CP_anti, stopping_time_CP_anti, dt_sde, r);

	for (size_t i = 0; i < simulated_price.nb_rows(); i++)
	{

		temp = simulated_price(i, 0) - simulated_price_CP(i, 0) + ExpPriceClsForm;
		temp_anti = simulated_price_anti(i, 0) - simulated_price_CP_anti(i, 0) + ExpPriceClsForm;

		average_price(i, 0) = (temp + temp_anti) * 0.5;

		temp = 0.;
		temp_anti = 0;

	}

	MC_price = average_price.mean();

	MC_variance = average_price.variance();

};
