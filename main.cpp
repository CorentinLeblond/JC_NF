#include "sde.hpp"
#include "MonteCarlo.h"
#include <iostream>
#include <fstream>
#include <chrono> 
#include <exception>

int main(int argc, char* argv[])
{	
	
	using clock = std::chrono::steady_clock;
/////////////////////////////////////////////////////////////////////////////////////
//INPUT FROM USER

	double startTime = 0.;
	double endTime = 1.;
	size_t nbsteps = 252;
	size_t Nb_Assets = 3;
	double rate = 0.05;
	double K = 100;
	
	std::vector<std::vector<double>> Spot_vector = {{100},{120},{80}};

	std::vector<std::vector<double>> Sigma_vector = {{0.25},{0.20},{0.15}};

	std::vector<std::vector<double>> Correl_mat = {{1,-0.2,0.4},
												   {-0.2,1,0.4},
												   {0.4,0.4,1}};
												   
	std::vector<std::vector<double>> Weights_mat ={{0.3,0.5,0.2}};
	
	matrix spot_m(Spot_vector);
	matrix Sigma(Sigma_vector);
	matrix Correl(Correl_mat);
	matrix Weights(Weights_mat);
	matrix CovarMatrix = VarCovarMatrix(Sigma, Correl);
	double df = exp(-rate * endTime);
	double dt = (endTime - startTime) / nbsteps;
/////////////////////////////////////////////////////////////////////////////////////
//Catching errors in inputs

// std::cout<<"nb rows "<< spot_m.nb_rows()<<std::endl;
	try
	{
		double sum_test = 0.;
		bool test_sym;
		
		//test Valeurs
		if(startTime<0.)
			throw std::exception("Start Time must be today or later");
		if(endTime<0. || endTime<startTime)
			throw std::exception("End Time must be later then start time and positive");
		if(nbsteps<=0.)
			throw std::exception("Number of steps must be strictly positive");
		if(Nb_Assets<=0.)
			throw std::exception("Number of assets must be strictly positive");
		
		for(size_t i=0;i<Sigma.nb_rows();++i)
		{
			if(Sigma(i,0)<0.)
				throw std::exception("Volatilities must be positive");			
		}
		for(size_t i=0;i<Weights.nb_cols();++i)
		{
			sum_test+=Weights(0,i);
		}		
		if(sum_test != 1.)
			throw std::exception("Sum of weights must be equal to 1");	
		
		for(size_t i=0;i<Correl.nb_rows();++i)
		{
			if(Correl(i,i)!= 1.)
				throw std::exception("diagonal elements of the correlation matrix must be 1s");			
		}
		for(size_t i=0;i<Correl.nb_rows();++i)
		{
			for(size_t j=0;j<Correl.nb_rows();++j)
			{
				if(Correl(i,j)< -1. || Correl(i,j)> 1.  )
					throw std::exception("All elements of the correlation matrix must be between -1 and 1");
			}			
		}
		
		//test Dimensions
		if(Nb_Assets!=spot_m.nb_rows())
			throw std::exception("The number of initial spot prices must be equal to the number of assets");
		if(Nb_Assets!=Sigma.nb_rows())
			throw std::exception("The number of volatilities must be equal to the number of assets");
		if(Nb_Assets!=Correl.nb_rows() && Nb_Assets!=Correl.nb_cols())
			throw std::exception("The correlation matrix must be a square matrix of dimension the number of assets");
		if(Nb_Assets!=Weights.nb_cols())
			throw std::exception("The number of weights must be equal to the number of assets");

		//test symmetry
		test_sym = testsymm(Correl);
		if(test_sym==false)
			throw std::exception("correlation matrix must be symmetric");
		
	}
	catch(std::exception const& e)
	{	
	   std::cout << "ERREUR : " << e.what() << std::endl;
	}	
	

/////////////////////////////////////////////////////////////////////////////////////
	UniformGenerator* ugen = new EcuyerCombined();
	
	Normal* ngen;
	
	try
	{
	ngen = new NormalBoxMuller(ugen, 0., 1.);
	}
	catch(std::exception const& e)
	{	
	   std::cout << "ERREUR : " << e.what() << std::endl;
	}
	
	GaussianVector* gvec;
	
	double determinant = determinantOfMatrix(CovarMatrix,CovarMatrix.nb_rows());
	
	if(determinant == 0.)
	{
		gvec = new GaussianVectorDiag(ngen,  Sigma,  Correl, CovarMatrix);
	}
	else
	{
		gvec = new GaussianVectorCholesky(ngen,  Sigma,  Correl, CovarMatrix); //More efficient method
	}
	
	RandomProcess* path = new BSEulerND(gvec,spot_m,rate);
	PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);

	//MONTE CARLO
	double tolerated_error = 0.01;
	size_t CompMC_Nb_Simulation = 2000;
	
	EuropeanBasket MC(CompMC_Nb_Simulation, bsktcall, path);
	
	//First Run: Companion MC to determinate the optimal number of simulation
	
	clock::time_point startComp = clock::now();
	MC.Simulate(startTime,endTime,nbsteps);
	clock::time_point endComp = clock::now();	
	clock::duration execution_timeComp = endComp - startComp;

	MC.OptimalNbSimul(tolerated_error);
	
	size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	std::cout << "EXECUTION TIME COMP: " << std::chrono::duration 
	<double,std::ratio<1>> (execution_timeComp).count() << std::endl;
	
	delete ugen;
	delete ngen;
	delete gvec;
	delete path;
	delete bsktcall;

	
//1. NO VARIANCE REDUCTION

	// RANDOM NUMBER GENERATION
	// UniformGenerator* ugen = new EcuyerCombined();
	// Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	// GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	// //RANDOM PROCESS
	// RandomProcess* path = new BSEulerND(corrG,spot_m,rate);

	// //PAYOFF
	// PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);

	// //MONTE CARLO
	// double tolerated_error = 0.01;
	// size_t CompMC_Nb_Simulation = 10000;
	
	// EuropeanBasket MC(CompMC_Nb_Simulation, bsktcall, path);
	
	// //First Run: Companion MC to determinate the optimal number of simulation
	
	// clock::time_point startComp = clock::now();
	// MC.Simulate(startTime,endTime,nbsteps);
	// clock::time_point endComp = clock::now();	
	// clock::duration execution_timeComp = endComp - startComp;

	// MC.OptimalNbSimul(tolerated_error);
	
	// size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	// std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	// std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	// std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	// std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	// std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	// std::cout << "EXECUTION TIME COMP: " << std::chrono::duration 
	// <double,std::ratio<1>> (execution_timeComp).count() << std::endl;


/////////////////////////////////////////////////////////////////////////////////////

//2.QUASI MC WITH SOBOL SEQUENCE
	// //RANDOM NUMBER GENERATION
	// UniformGenerator* ugen = new Sobol();
	// Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	// GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	// //RANDOM PROCESS
	// RandomProcess* path = new BSEulerND(corrG,spot_m,rate);

	// //PAYOFF
	// PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);

	// //MONTE CARLO
	// double tolerated_error = 0.01;
	// size_t CompMC_Nb_Simulation = 10000;
	
	// EuropeanBasket MC(CompMC_Nb_Simulation, bsktcall, path);
	
	// //First Run: Companion MC to determinate the optimal number of simulation
	
	// clock::time_point startComp = clock::now();
	// MC.Simulate(startTime,endTime,nbsteps);
	// clock::time_point endComp = clock::now();	
	// clock::duration execution_timeComp = endComp - startComp;

	// MC.OptimalNbSimul(tolerated_error);
	
	// size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	// std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	// std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	// std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	// std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	// std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	// std::cout << "EXECUTION TIME COMP: " << std::chrono::duration 
	// <double,std::ratio<1>> (execution_timeComp).count() << std::endl;



/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

//2.5.QUASI MC WITH VDC SEQUENCE
	//RANDOM NUMBER GENERATION
	// UniformGenerator* ugen = new VanDerCorput(2,1);
	// UniformGenerator* ugen2 = new VanDerCorput(3,1);
	// Normal* ngen = new NormalBoxMullerVDC(ugen,ugen2, 0., 1.);
	// GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	// //RANDOM PROCESS
	// RandomProcess* path = new BSEulerND(corrG,spot_m,rate);

	// //PAYOFF
	// PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);

	// //MONTE CARLO
	// double tolerated_error = 0.01;
	// size_t CompMC_Nb_Simulation = 3000;
	
	// EuropeanBasket MC(CompMC_Nb_Simulation, bsktcall, path);
	
	// //First Run: Companion MC to determinate the optimal number of simulation
	
	// clock::time_point startComp = clock::now();
	// MC.Simulate(startTime,endTime,nbsteps);
	// clock::time_point endComp = clock::now();	
	// clock::duration execution_timeComp = endComp - startComp;

	// MC.OptimalNbSimul(tolerated_error);
	
	// size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	// std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	// std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	// std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	// std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	// std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	// std::cout << "EXECUTION TIME COMP: " << std::chrono::duration 
	// <double,std::ratio<1>> (execution_timeComp).count() << std::endl;

// delete ugen2;

/////////////////////////////////////////////////////////////////////////////////////
//3.ANTHITETIC MC
	//RANDOM NUMBER GENERATION
	// UniformGenerator* ugen = new EcuyerCombined();
	// Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	// GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	// //RANDOM PROCESS
	// RandomProcess* path = new BSEulerNDAntithetic(corrG,spot_m,rate);

	// //PAYOFF
	// PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);

	// //MONTE CARLO
	// double tolerated_error = 0.01;
	// size_t CompMC_Nb_Simulation = 3000;
	
	// EuropeanBasket_Antithetic MC(CompMC_Nb_Simulation, bsktcall, path);
	
	// //First Run: Companion MC to determinate the optimal number of simulation
	
	// clock::time_point startComp = clock::now();
	// MC.Simulate(startTime,endTime,nbsteps);
	// clock::time_point endComp = clock::now();
	// clock::duration execution_timeComp = endComp - startComp;
	
	// MC.OptimalNbSimul(tolerated_error);
	
	// size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	// std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	// std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	// std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	// std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	// std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	// std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;


/////////////////////////////////////////////////////////////////////////////////////

//4.CV MC
	//RANDOM NUMBER GENERATION
	
	// UniformGenerator* ugen = new EcuyerCombined();
	
	// Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	// GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	// // UniformGenerator* ugen = new Sobol();
	
	// //RANDOM PROCESS
	// RandomProcess* path = new BSEulerND(corrG,spot_m,rate);

	// //PAYOFF
	// PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);
	
	// PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(Weights, spot_m,K);
	// ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall
	// (Weights, spot_m, CovarMatrix, K,rate, endTime);
	// double CV_closedprice = CFbaskt->operator()(spot_m,df);
	
	// //MONTE CARLO
	// double tolerated_error = 0.01;
	// size_t CompMC_Nb_Simulation = 2000;
	
	// EuropeanBasket_controlvariable MC(CompMC_Nb_Simulation,
	// bsktcall,bsktcallCV, path,CV_closedprice);
	
	// //First Run: Companion MC to determinate the optimal number of simulation
	
	// clock::time_point startComp = clock::now();
	// MC.Simulate(startTime,endTime,nbsteps);
	// clock::time_point endComp = clock::now();
	// clock::duration execution_timeComp = endComp - startComp;
	
	// //dynamics.Simulate(startTime, endTime, nbsteps);
	// //matrix chemin = dynamics.GetAllPaths();
	// //chemin.Print();

	// UniformGenerator* ugen2 = new EcuyerCombined();
	// RandomGenerator* ngen2 = new NormalBoxMuller(ugen2, 0., 1.);
	// // RandomProcess* chemin = new BSEuler1D(ngen2, spot, rate, vol);

	
	// MC.OptimalNbSimul(tolerated_error);
	
	// size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	// std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	// std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	// std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	// std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	// std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	// std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;
	// delete CFbaskt;
	// delete bsktcallCV;

/////////////////////////////////////////////////////////////////////////////////////



//5.ANTHITETIC + SOBOL MC
	//RANDOM NUMBER GENERATION
	// UniformGenerator* ugen = new Sobol();
	// Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	// GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	// //RANDOM PROCESS
	// RandomProcess* path = new BSEulerNDAntithetic(corrG,spot_m,rate);

	// //PAYOFF
	// PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);

	// //MONTE CARLO
	// double tolerated_error = 0.01;
	// size_t CompMC_Nb_Simulation = 2000;
	
	// EuropeanBasket_Antithetic MC(CompMC_Nb_Simulation, bsktcall, path);
	
	// //First Run: Companion MC to determinate the optimal number of simulation
	
	// clock::time_point startComp = clock::now();
	// MC.Simulate(startTime,endTime,nbsteps);
	// clock::time_point endComp = clock::now();
	// clock::duration execution_timeComp = endComp - startComp;
	
	
	// MC.OptimalNbSimul(tolerated_error);
	
	// size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	// std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	// std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	// std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	// std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	// std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	// std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;


/////////////////////////////////////////////////////////////////////////////////////

//6. QUASI + CV MC
	// RANDOM NUMBER GENERATION
	// UniformGenerator* ugen = new Sobol();
	// Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	// GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	// // RANDOM PROCESS
	// RandomProcess* path = new BSEulerND(corrG,spot_m,rate);

	// // PAYOFF
	// PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);
	// PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(Weights, spot_m,K);
	// ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(Weights, spot_m, CovarMatrix, K,rate, endTime);

	// double CV_closedprice = CFbaskt->operator()(spot_m,df);
	// // std::cout << "closed formula for the basket option" << CV_closedprice << std::endl;
	// // MONTE CARLO
	// double tolerated_error = 0.01;
	// size_t CompMC_Nb_Simulation = 2000;
	
	// EuropeanBasket_controlvariable MC(CompMC_Nb_Simulation, bsktcall,bsktcallCV, path,CV_closedprice);
	
	// // First Run: Companion MC to determinate the optimal number of simulation
	
	// clock::time_point startComp = clock::now();
	// MC.Simulate(startTime,endTime,nbsteps);
	// clock::time_point endComp = clock::now();
	// clock::duration execution_timeComp = endComp - startComp;
	
	
	// MC.OptimalNbSimul(tolerated_error);
	
	// size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	// std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	// std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	// std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	// std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	// std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	// std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;
	// delete CFbaskt;
	// delete bsktcallCV;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//7. Anti + CV MC
	//RANDOM NUMBER GENERATION
	// UniformGenerator* ugen = new EcuyerCombined();
	// Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	// GaussianVectorCholesky* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	// //RANDOM PROCESS
	// BSEulerNDAntithetic* path = new BSEulerNDAntithetic(corrG,spot_m,rate);

	// //PAYOFF
	// PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);
	// PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(Weights, spot_m,K);
	// ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(Weights, spot_m, CovarMatrix, K,rate, endTime);

	// double CV_closedprice = CFbaskt->operator()(spot_m,df);
	// // std::cout << "closed formula for the basket option" << CV_closedprice << std::endl;
	// //MONTE CARLO
	// double tolerated_error = 0.01;
	// size_t CompMC_Nb_Simulation = 3000;
	
	// EuropeanBasket_Antithetic_CV MC(CompMC_Nb_Simulation,
	// bsktcall,bsktcallCV, path,CV_closedprice);
	
	// //First Run: Companion MC to determinate the optimal number of simulation
	
	// clock::time_point startComp = clock::now();
	// MC.Simulate(startTime,endTime,nbsteps);
	// clock::time_point endComp = clock::now();
	// clock::duration execution_timeComp = endComp - startComp;
	
	
	// MC.OptimalNbSimul(tolerated_error);
	
	// size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	// std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	// std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	// std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	// std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	// std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	// std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;
	// delete CFbaskt;
	// delete bsktcallCV;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//8. 3 TECHNIQUES USED
	//RANDOM NUMBER GENERATION
	// UniformGenerator* ugen = new Sobol();
	// Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	// GaussianVectorCholesky* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	// //RANDOM PROCESS
	// BSEulerNDAntithetic* path = new BSEulerNDAntithetic(corrG,spot_m,rate);

	// //PAYOFF
	// PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);
	// PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(Weights, spot_m,K);
	// ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(Weights, spot_m, CovarMatrix, K,rate, endTime);

	// double CV_closedprice = CFbaskt->operator()(spot_m,df);
	// // std::cout << "closed formula for the basket option" << CV_closedprice << std::endl;
	// //MONTE CARLO
	// double tolerated_error = 0.01;
	// size_t CompMC_Nb_Simulation = 2000;
	
	// EuropeanBasket_Antithetic_CV MC(CompMC_Nb_Simulation, bsktcall,bsktcallCV, path,CV_closedprice);
	
	// //First Run: Companion MC to determinate the optimal number of simulation
	
	// clock::time_point startComp = clock::now();
	// MC.Simulate(startTime,endTime,nbsteps);
	// clock::time_point endComp = clock::now();
	// clock::duration execution_timeComp = endComp - startComp;
	
	
	// MC.OptimalNbSimul(tolerated_error);
	
	// size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	// std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	// std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	// std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	// std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	// std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	// std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;
	// delete CFbaskt;
	// delete bsktcallCV;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*	
	std::cout << "PRICING"<< std::endl;
	clock::time_point start = clock::now();
	
	//Second Run: MC to compute to price of the option with required confidence level
	MC.Simulate(startTime,endTime,Optimal_Nb_Simulation);
	
	clock::time_point end = clock::now(); 
	
	clock::duration execution_time = end - start;
	
	std::cout << "EXECUTION TIME: " << std::chrono::duration <double,std::ratio<1>> (execution_time).count() << std::endl;
	
	double priceMC = MC.GetPrice(rate, endTime);
	double varMC = MC.GetVariance();

	std::cout << "OPTION PRICE: " << priceMC << std::endl;
	std::cout << "MC VARIANCE: " << varMC << std::endl;
	std::cout << "NB SIMULATION: " << Optimal_Nb_Simulation << std::endl;
*/

	// delete ugen;
	// delete ngen;
	// delete corrG;
	// delete path;
	// delete bsktcall;

	// EuropeanBasket_controlvariable CVMC(N, bsktcall, bsktcallCV, path_cv);
	// CVMC.Simulate(startTime, endTime, nbsteps);
	// double priceCVMC = CVMC.GetPrice(rate, endTime);

	// double varVMC = VMC.GetVariance();
	// double varMC = MC.GetVariance();
	// double varCVMC = CVMC.GetVariance();

	// std::cout << "price Vanilla 1D MC " << priceVMC << std::endl;
	// std::cout << "price Vanilla CV MC " << priceCVMC << std::endl;
	// std::cout << "price Vanilla MC for Basket " << priceMC << std::endl;

	// std::cout << "variance Vanilla 1D MC " << varVMC << std::endl;
	// std::cout << "variance Vanilla CV MC " << varCVMC << std::endl;
	// std::cout << "variance Vanilla MC for Basket " << varMC << std::endl;
	
	// ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(W, spot_m, CovarMatrix, K,
		// rate, endTime);

	//std::vector<std::vector<double>> Spot_vector_maturity = { {130},{110},{91} };
	//matrix endspot(Spot_vector_maturity);
	// double df = exp(-rate * endTime);
	// std::cout << "closed formula for the basket option" << CFbaskt->operator()(spot_m,df) << std::endl;
	
	// std::cout << "debut print vector" << std::endl;
	// for(size_t i = 0; i< nbsteps;++i)
	// {
		// std::cout << "Step " << i << ", Spot= " << path1->GetState(i*dt) << ",";
	// }
	// std::cout << std::endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Test GaussianVector generation

	//matrix Sigma(Sigma_vector);
	//matrix Nu(Mu_vector);
	//matrix Correl(Correl_mat);

	//matrix CovarMatrix = VarCovarMatrix(Sigma,Correl);
	//
	//Normal* ngnr = new NormalBoxMuller(ptr,0.,1.);
	//
	//// bool testinvertible = isInvertible(CovarMatrix,CovarMatrix.nb_rows());
	//// std::cout<< "test invertivle returns: " << testinvertible << std::endl;
	//GaussianVector* gvec;
	//
	//// if(testinvertible == 1)
	//gvec = new GaussianVectorCholesky(ngnr,Nu,Sigma,Correl,CovarMatrix);
	//// else
	//	// gvec = new GaussianVectorDiag(ngnr,Nu,Sigma,Correl,CovarMatrix);

	//matrix output = gvec->CorrelatedGaussianVector();
	//std::cout << "Output" <<std::endl;
	//output.Print();
	//matrix output2 = gvec->CorrelatedGaussianVector();
	//std::cout << "Output2" <<std::endl;
	//output2.Print();

//////////////////////////// US Part ////////////////////////////////////////////////////

std::vector<basis_functions*> basefunc_Laguerre;

basefunc_Laguerre.push_back(new Poly_Laguerre(0));
basefunc_Laguerre.push_back(new Poly_Laguerre(1));
basefunc_Laguerre.push_back(new Poly_Laguerre(2));


UniformGenerator* ugen = new EcuyerCombined();
Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);

GaussianVector* vectorG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican = new BSEulerND(vectorG, spot_m, rate);

//MONTE CARLO
AmericanMonteCarlo_basket USMC(N, bsktcall, BSamerican, basefunc_Laguerre);

clock::time_point start_US = clock::now();
USMC.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US = clock::now();
clock::duration execution_timeUS = end_US - start_US;
double price_US = USMC.GetPrice();
double var_US = USMC.GetVariance();

USMC.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation = USMC.GetNbSimul();

std::cout << " US Case - exec time for vanilla MonteCarlo : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS).count() << std::endl;
std::cout << " US Case - price for vanilla MonteCarlo : " << price_US << std::endl;
std::cout << " US Case - variance for vanilla MonteCarlo : " << var_US << std::endl;
std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo : " << Optimal_Nb_Simulation << std::endl;

delete ugen;
delete ngen;
delete vectorG;
delete BSamerican;

//// need to choose the polynome basis with which one will approxime the continuation value of the option
//	std::vector<basis_functions*> basefunc_Hermite;
//
//	basefunc_Hermite.push_back(new Poly_Hermite(0));
//	basefunc_Hermite.push_back(new Poly_Hermite(1));
//	basefunc_Hermite.push_back(new Poly_Hermite(2));
//
//	std::vector<basis_functions*> basefunc_Laguerre;
//
//	basefunc_Laguerre.push_back(new Poly_Laguerre(0));
//	basefunc_Laguerre.push_back(new Poly_Laguerre(1));
//	basefunc_Laguerre.push_back(new Poly_Laguerre(2));
//
//	std::vector<basis_functions*> basefunc_Simple;
//
//	basefunc_Simple.push_back(new polynome_simple(0));
//	basefunc_Simple.push_back(new polynome_simple(1));
//	basefunc_Simple.push_back(new polynome_simple(2));
//
//	//Create the payoff of the options 
//
//	PayOffBasket* bsktcall = new PayOffBasketCall(W, spot_m, 100.);
//
//	/////////////////// NO VARIANCE REDUCTION /////////////////////////////
//
//	////////////////// 1.1) Vanille (Ecuyer Combined - Laguerre Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* ugen = new EcuyerCombined();
//	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
//
//	GaussianVector* vectorG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican = new BSEulerND(vectorG, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC(N, bsktcall, BSamerican, basefunc_Laguerre);
//
//	clock::time_point start_US = clock::now();
//	USMC.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US = clock::now(); 
//	clock::duration execution_timeUS = end_US - start_US; 
//	double price_US = USMC.GetPrice();
//	double var_US = USMC.GetVariance();
//
//	USMC.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation = USMC.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo : " << price_US << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo : " << var_US << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo : " << Optimal_Nb_Simulation << std::endl;
//
//	delete ugen;
//	delete ngen;
//	delete vectorG;
//	delete BSamerican;
//
//	////////////////// 1.1.2) Vanille (VDC - Laguerre Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* vdc_gen = new VanDerCorput(2,1);
//	UniformGenerator* vdc_gen_alias = new VanDerCorput(11, 1);
//	Normal* ngen_vdc = new NormalBoxMullerVDC(vdc_gen, vdc_gen_alias, 0., 1.);
//
//	GaussianVector* vectorG_vdc = new GaussianVectorCholesky(ngen_vdc, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc = new BSEulerND(vectorG_vdc, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_vdc(N, bsktcall, BSamerican_vdc, basefunc_Laguerre);
//
//	clock::time_point start_US_vdc = clock::now();
//	USMC_vdc.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc = clock::now();
//	clock::duration execution_timeUS_vdc = end_US_vdc - start_US_vdc;
//	double price_US_vdc = USMC_vdc.GetPrice();
//	double var_US_vdc = USMC_vdc.GetVariance();
//
//	USMC_vdc.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc = USMC_vdc.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with VDC sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with VDC sequence : " << price_US_vdc << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo  with VDC sequence : " << var_US_vdc << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo  with VDC sequence : " << Optimal_Nb_Simulation_vdc << std::endl;
//
//	delete vdc_gen;
//	delete vdc_gen_alias;
//	delete ngen_vdc;
//	delete vectorG_vdc;
//	delete BSamerican_vdc;
//
//	////////////////// 1.1.3) Vanille (Sobol - Laguerre Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* sobol_ugen = new Sobol();
//	Normal* ngen_sobol = new NormalBoxMuller(sobol_ugen, 0., 1.);
//
//	GaussianVector* vectorG_sobol = new GaussianVectorCholesky(ngen_sobol, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol = new BSEulerND(vectorG_sobol, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_sobol(N, bsktcall, BSamerican_sobol, basefunc_Laguerre);
//
//	clock::time_point start_US_sobol = clock::now();
//	USMC_sobol.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol = clock::now();
//	clock::duration execution_timeUS_sobol = end_US_sobol - start_US_sobol;
//	double price_US_sobol = USMC_sobol.GetPrice();
//	double var_US_sobol = USMC_sobol.GetVariance();
//
//	USMC_sobol.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol = USMC_sobol.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with Sobol sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with Sobol sequence : " << price_US_sobol << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo  with Sobol sequence : " << var_US_sobol << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo  with Sobol sequence : " << Optimal_Nb_Simulation_sobol << std::endl;
//
//	delete sobol_ugen;
//	delete ngen_sobol;
//	delete vectorG_sobol;
//	delete BSamerican_sobol;
//
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	////////////////// 1.2) Vanille (Ecuyer Combined - Hermite Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* ugen_H = new EcuyerCombined();
//	Normal* ngen_H = new NormalBoxMuller(ugen_H, 0., 1.);
//
//	GaussianVector* vectorG_H = new GaussianVectorCholesky(ngen_H, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_H = new BSEulerND(vectorG_H, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_Hermite(N, bsktcall, BSamerican_H, basefunc_Hermite);
//
//	clock::time_point start_US_H = clock::now();
//	USMC_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_H = clock::now();
//	clock::duration execution_timeUS_H = end_US_H - start_US_H;
//	double price_US_H = USMC_Hermite.GetPrice();
//	double var_US_H = USMC_Hermite.GetVariance();
//
//	USMC_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_H = USMC_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with Hermite polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_H).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with Hermite polynome : " << price_US_H << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo with Hermite polynome : " << var_US_H << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo with Hermite polynome : " << Optimal_Nb_Simulation_H << std::endl;
//
//	delete ugen_H;
//	delete ngen_H;
//	delete vectorG_H;
//	delete BSamerican_H;
//
//	////////////////// 1.2.2) Vanille (VDC - Hermite Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* vdc_gen_H = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_H = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_H = new NormalBoxMullerVDC(vdc_gen_H, vdc_gen_alias_H, 0., 1.);
//
//	GaussianVector* vectorG_vdc_H = new GaussianVectorCholesky(ngen_vdc_H, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_H = new BSEulerND(vectorG_vdc_H, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_vdc_Hermite(N, bsktcall, BSamerican_vdc_H, basefunc_Hermite);
//
//	clock::time_point start_US_vdc_H = clock::now();
//	USMC_vdc_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_H = clock::now();
//	clock::duration execution_timeUS_vdc_H = end_US_vdc_H - start_US_vdc_H;
//	double price_US_vdc_H = USMC_vdc_Hermite.GetPrice();
//	double var_US_vdc_H = USMC_vdc_Hermite.GetVariance();
//
//	USMC_vdc_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_Hermite = USMC_vdc_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with VDC sequence and Hermite Polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_H).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with VDC sequence and Hermite Polynomes : " << price_US_vdc_H << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo  with VDC sequence and Hermite Polynomes : " << var_US_vdc_H << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo  with VDC sequence and Hermite Polynomes : " << Optimal_Nb_Simulation_vdc_Hermite << std::endl;
//
//	delete vdc_gen_H;
//	delete vdc_gen_alias_H;
//	delete ngen_vdc_H;
//	delete vectorG_vdc_H;
//	delete BSamerican_vdc_H;
//
//	////////////////// 1.2.3) Vanille (Sobol - Hermite Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* sobol_ugen_H = new Sobol();
//	Normal* ngen_sobol_H = new NormalBoxMuller(sobol_ugen_H, 0., 1.);
//
//	GaussianVector* vectorG_sobol_H = new GaussianVectorCholesky(ngen_sobol_H, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_H = new BSEulerND(vectorG_sobol_H, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_sobol_Hermite(N, bsktcall, BSamerican_sobol_H, basefunc_Hermite);
//
//	clock::time_point start_US_sobol_H = clock::now();
//	USMC_sobol_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_H = clock::now();
//	clock::duration execution_timeUS_sobol_H = end_US_sobol_H - start_US_sobol_H;
//	double price_US_sobol_H = USMC_sobol_Hermite.GetPrice();
//	double var_US_sobol_H = USMC_sobol_Hermite.GetVariance();
//
//	USMC_sobol_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_H = USMC_sobol_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with Sobol sequence and Hermite Polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_H).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with Sobol sequence and Hermite Polynomes : " << price_US_sobol_H << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo  with Sobol sequence and Hermite Polynomes : " << var_US_sobol_H << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo  with Sobol sequence and Hermite Polynomes : " << Optimal_Nb_Simulation_sobol_H << std::endl;
//
//	delete sobol_ugen_H;
//	delete ngen_sobol_H;
//	delete vectorG_sobol_H;
//	delete BSamerican_sobol_H;
//
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	////////////////// 1.3) Vanille (Ecuyer Combined - Simples Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* ugen_S = new EcuyerCombined();
//	Normal* ngen_S = new NormalBoxMuller(ugen_S, 0., 1.);
//
//	GaussianVector* vectorG_S = new GaussianVectorCholesky(ngen_S, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_S = new BSEulerND(vectorG_S, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_Simples(N, bsktcall, BSamerican_S, basefunc_Simple);
//
//	clock::time_point start_US_S = clock::now();
//	USMC_Simples.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_S = clock::now();
//	clock::duration execution_timeUS_S = end_US_S - start_US_S;
//	double price_US_S = USMC_Simples.GetPrice();
//	double var_US_S = USMC_Simples.GetVariance();
//
//	USMC_Simples.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_S = USMC_Simples.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with simples polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_S).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with simples polynome : " << price_US_S << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo with simples polynome : " << var_US_S << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo with simples polynome : " << Optimal_Nb_Simulation_S << std::endl;
//
//	delete ugen_S;
//	delete ngen_S;
//	delete vectorG_S;
//	delete BSamerican_S;
//
//	////////////////// 1.3.2) Vanille (VDC - simples Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* vdc_gen_S = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_S = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_S = new NormalBoxMullerVDC(vdc_gen_S, vdc_gen_alias_S, 0., 1.);
//
//	GaussianVector* vectorG_vdc_S = new GaussianVectorCholesky(ngen_vdc_S, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_S = new BSEulerND(vectorG_vdc_S, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_vdc_Simples(N, bsktcall, BSamerican_vdc_S, basefunc_Simple);
//
//	clock::time_point start_US_vdc_S = clock::now();
//	USMC_vdc_Simples.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_S = clock::now();
//	clock::duration execution_timeUS_vdc_S = end_US_vdc_S - start_US_vdc_S;
//	double price_US_vdc_S = USMC_vdc_Simples.GetPrice();
//	double var_US_vdc_S = USMC_vdc_Simples.GetVariance();
//
//	USMC_vdc_Simples.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_Simples = USMC_vdc_Simples.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with VDC sequence and simples Polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_S).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with VDC sequence and simples Polynomes : " << price_US_vdc_S << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo  with VDC sequence and simples Polynomes : " << var_US_vdc_S << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo  with VDC sequence and simples Polynomes : " << Optimal_Nb_Simulation_vdc_Simples << std::endl;
//
//	delete vdc_gen_S;
//	delete vdc_gen_alias_S;
//	delete ngen_vdc_S;
//	delete vectorG_vdc_S;
//	delete BSamerican_vdc_S;
//
//	////////////////// 1.3.3) Vanille (Sobol - simples Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* sobol_ugen_S = new Sobol();
//	Normal* ngen_sobol_S = new NormalBoxMuller(sobol_ugen_S, 0., 1.);
//
//	GaussianVector* vectorG_sobol_S = new GaussianVectorCholesky(ngen_sobol_S, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_S = new BSEulerND(vectorG_sobol_S, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_sobol_Simples(N, bsktcall, BSamerican_sobol_S, basefunc_Simple);
//
//	clock::time_point start_US_sobol_S = clock::now();
//	USMC_sobol_Simples.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_S = clock::now();
//	clock::duration execution_timeUS_sobol_S = end_US_sobol_S - start_US_sobol_S;
//	double price_US_sobol_S = USMC_sobol_Simples.GetPrice();
//	double var_US_sobol_S = USMC_sobol_Simples.GetVariance();
//
//	USMC_sobol_Simples.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_S = USMC_sobol_Simples.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with Sobol sequence and simples Polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_S).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with Sobol sequence and simples Polynomes : " << price_US_sobol_S << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo  with Sobol sequence and simples Polynomes : " << var_US_sobol_S << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo  with Sobol sequence and simples Polynomes : " << Optimal_Nb_Simulation_sobol_S << std::endl;
//
//	delete sobol_ugen_S;
//	delete ngen_sobol_S;
//	delete vectorG_sobol_S;
//	delete BSamerican_sobol_S;
//
/////////////////////////////////////// VARIANCE REDUCTION TECHNIQUES ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//	//////////////////////////// CONTROL VARIABLE /////////////////////////////////////////////////////////////////////
//
//
//	///create both the control variate payoff and its close_formula price in order to proceed to the contrl variate technique
//
//	PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(W, spot_m, 100.);
//
//	ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(W, spot_m, CovarMatrix, K,
//		rate, endTime);
//
//	double df = exp(-rate * endTime);
//	double price_cf = CFbaskt->operator()(spot_m, df);
//
////////////////////////// 2.1) Control Variable (Ecuyer Combined - Laguerre polynome n = 2) ////////////////////////////
//
//
//	UniformGenerator* ugen_VC = new EcuyerCombined();
//	Normal* ngen_VC = new NormalBoxMuller(ugen_VC, 0., 1.);
//
//	GaussianVector* vectorG_VC = new GaussianVectorCholesky(ngen_VC, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_VC = new BSEulerND(vectorG_VC, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_VC(N, bsktcall, bsktcallCV, BSamerican_VC,
//		basefunc_Laguerre, price_cf);
//
//	clock::time_point start_US_VC = clock::now();
//	USMC_VC.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_VC = clock::now();
//	clock::duration execution_timeUS_VC = end_US_VC - start_US_VC;
//	double price_US_VC = USMC_VC.GetPrice();
//	double var_US_VC = USMC_VC.GetVariance();
//
//	USMC_VC.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_VC = USMC_VC.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_VC).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo : " << price_US_VC << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo : " << var_US_VC << std::endl;
//	std::cout << " US Case - optimal number of simulations for control variable MonteCarlo : " << Optimal_Nb_Simulation_VC << std::endl;
//
//	delete ugen_VC;
//	delete ngen_VC;
//	delete vectorG_VC;
//	delete BSamerican_VC;
//
//	////////////////// 2.1.2) Control variable (VDC - Laguerre Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* vdc_gen_VC = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_VC = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_VC = new NormalBoxMullerVDC(vdc_gen_VC, vdc_gen_alias_VC, 0., 1.);
//
//	GaussianVector* vectorG_vdc_VC = new GaussianVectorCholesky(ngen_vdc_VC, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_VC = new BSEulerND(vectorG_vdc_VC, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_VC_vdc(N, bsktcall, bsktcallCV, BSamerican_vdc_VC,
//		basefunc_Laguerre, price_cf);
//
//	clock::time_point start_US_vdc_VC = clock::now();
//	USMC_VC_vdc.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_VC = clock::now();
//	clock::duration execution_timeUS_vdc_VC = end_US_vdc_VC - start_US_vdc_VC;
//	double price_US_vdc_VC = USMC_VC_vdc.GetPrice();
//	double var_US_vdc_VC = USMC_VC_vdc.GetVariance();
//
//	USMC_VC_vdc.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_VC = USMC_VC_vdc.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo with VDC sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_VC).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo with VDC sequence : " << price_US_vdc_VC << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo  with VDC sequence : " << var_US_vdc_VC << std::endl;
//	std::cout << " US Case - optimal number of simulations for control variable MonteCarlo with VDC sequence : " << Optimal_Nb_Simulation_vdc_VC << std::endl;
//
//	delete vdc_gen_VC;
//	delete vdc_gen_alias_VC;
//	delete ngen_vdc_VC;
//	delete vectorG_vdc_VC;
//	delete BSamerican_vdc_VC;
//
//	////////////////// 2.1.3) Control variable (Sobol - Laguerre Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* sobol_ugen_CV = new Sobol();
//	Normal* ngen_sobol_CV = new NormalBoxMuller(sobol_ugen_CV, 0., 1.);
//
//	GaussianVector* vectorG_sobol_CV = new GaussianVectorCholesky(ngen_sobol_CV, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_CV = new BSEulerND(vectorG_sobol_CV, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_sobol_CV(N, bsktcall, bsktcallCV, BSamerican_sobol_CV, basefunc_Laguerre, price_cf);
//
//	clock::time_point start_US_sobol_CV = clock::now();
//	USMC_sobol_CV.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_CV = clock::now();
//	clock::duration execution_timeUS_sobol_CV = end_US_sobol_CV - start_US_sobol_CV;
//	double price_US_sobol_CV = USMC_sobol_CV.GetPrice();
//	double var_US_sobol_CV = USMC_sobol_CV.GetVariance();
//
//	USMC_sobol_CV.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_CV = USMC_sobol_CV.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo with Sobol sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_CV).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo with Sobol sequence : " << price_US_sobol_CV << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo with Sobol sequence : " << var_US_sobol_CV << std::endl;
//	std::cout << " US Case - optimal number of simulations forcontrol variable  MonteCarlo  with Sobol sequence : " << Optimal_Nb_Simulation_sobol_CV << std::endl;
//
//	delete sobol_ugen_CV;
//	delete ngen_sobol_CV;
//	delete vectorG_sobol_CV;
//	delete BSamerican_sobol_CV;
//
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	////////////////// 2.2) Control variable (Ecuyer Combined - Hermite Polynomes n = 2) //////////////////////////////////////
//
//
//	UniformGenerator* ugen_VCH = new EcuyerCombined();
//	Normal* ngen_VCH = new NormalBoxMuller(ugen_VCH, 0., 1.);
//
//	GaussianVector* vectorG_VCH = new GaussianVectorCholesky(ngen_VCH, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_VCH = new BSEulerND(vectorG_VCH, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_VCH(N, bsktcall, bsktcallCV, BSamerican_VCH,
//		basefunc_Hermite, price_cf);
//
//	clock::time_point start_US_VCH = clock::now();
//	USMC_VCH.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_VCH = clock::now();
//	clock::duration execution_timeUS_VCH = end_US_VCH - start_US_VCH;
//	double price_US_VCH = USMC_VCH.GetPrice();
//	double var_US_VCH = USMC_VCH.GetVariance();
//
//	USMC_VCH.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_VCH = USMC_VCH.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo and Hermite polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_VCH).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo and Hermite polynomes  : " << price_US_VCH << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo and Hermite polynomes  : " << var_US_VCH << std::endl;
//	std::cout << " US Case - optimal number of simulations for control variable MonteCarlo and Hermite polynomes  : " << Optimal_Nb_Simulation_VCH << std::endl;
//
//	delete ugen_VCH;
//	delete ngen_VCH;
//	delete vectorG_VCH;
//	delete BSamerican_VCH;
//
//	////////////////// 2.2.2) Control variable (VDC - Hermite Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* vdc_gen_VCH = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_VCH = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_VCH = new NormalBoxMullerVDC(vdc_gen_VCH, vdc_gen_alias_VCH, 0., 1.);
//
//	GaussianVector* vectorG_vdc_VCH = new GaussianVectorCholesky(ngen_vdc_VCH, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_VCH = new BSEulerND(vectorG_vdc_VCH, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_VC_vdc_Hermite(N, bsktcall, bsktcallCV, BSamerican_vdc_VCH,
//		basefunc_Hermite, price_cf);
//
//	clock::time_point start_US_vdc_VCH = clock::now();
//	USMC_VC_vdc_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_VCH = clock::now();
//	clock::duration execution_timeUS_vdc_VCH = end_US_vdc_VCH - start_US_vdc_VCH;
//	double price_US_vdc_VCH = USMC_VC_vdc_Hermite.GetPrice();
//	double var_US_vdc_VCH = USMC_VC_vdc_Hermite.GetVariance();
//
//	USMC_VC_vdc_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_VCH = USMC_VC_vdc_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo with VDC sequence and Hermite polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_VCH).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo with VDC sequence and Hermite polynomes  : " << price_US_vdc_VCH << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo  with VDC sequence and Hermite polynomes : " << var_US_vdc_VCH << std::endl;
//	std::cout << " US Case - optimal number of simulations for control variable MonteCarlo with VDC sequence and Hermite polynomes : " << Optimal_Nb_Simulation_vdc_VCH << std::endl;
//
//	delete vdc_gen_VCH;
//	delete vdc_gen_alias_VCH;
//	delete ngen_vdc_VCH;
//	delete vectorG_vdc_VCH;
//	delete BSamerican_vdc_VCH;
//
//	////////////////// 2.2.3) Control variable (Sobol -Hermite Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* sobol_ugen_CVH = new Sobol();
//	Normal* ngen_sobol_CVH = new NormalBoxMuller(sobol_ugen_CVH, 0., 1.);
//
//	GaussianVector* vectorG_sobol_CVH = new GaussianVectorCholesky(ngen_sobol_CVH, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_CVH = new BSEulerND(vectorG_sobol_CVH, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_sobol_CV_Hermite(N, bsktcall, bsktcallCV,BSamerican_sobol_CVH, basefunc_Hermite, price_cf);
//
//	clock::time_point start_US_sobol_CVH = clock::now();
//	USMC_sobol_CV_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_CVH = clock::now();
//	clock::duration execution_timeUS_sobol_CVH = end_US_sobol_CVH - start_US_sobol_CVH;
//	double price_US_sobol_CVH = USMC_sobol_CV_Hermite.GetPrice();
//	double var_US_sobol_CVH = USMC_sobol_CV_Hermite.GetVariance();
//
//	USMC_sobol_CV_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_CVH = USMC_sobol_CV_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo with Sobol sequence and Hermite polynomes  : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_CVH).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo with Sobol sequence and Hermite polynomes  : " << price_US_sobol_CVH << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo with Sobol sequence and Hermite polynomes  : " << var_US_sobol_CVH << std::endl;
//	std::cout << " US Case - optimal number of simulations forcontrol variable  MonteCarlo  with Sobol sequence and Hermite polynomes  : " << Optimal_Nb_Simulation_sobol_CVH << std::endl;
//
//	delete sobol_ugen_CVH;
//	delete ngen_sobol_CVH;
//	delete vectorG_sobol_CVH;
//	delete BSamerican_sobol_CVH;
//
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	////////////////// 2.3) Controle variable (Ecuyer Combined - Simples Polynomes n = 2) //////////////////////////////////////
//
//
//	UniformGenerator* ugen_VCS = new EcuyerCombined();
//	Normal* ngen_VCS = new NormalBoxMuller(ugen_VCS, 0., 1.);
//
//	GaussianVector* vectorG_VCS = new GaussianVectorCholesky(ngen_VCS, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_VCS = new BSEulerND(vectorG_VCS, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_VCS(N, bsktcall, bsktcallCV, BSamerican_VCS,
//		basefunc_Simple, price_cf);
//
//	clock::time_point start_US_VCS = clock::now();
//	USMC_VCS.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_VCS = clock::now();
//	clock::duration execution_timeUS_VCS = end_US_VCS - start_US_VCS;
//	double price_US_VCS = USMC_VCS.GetPrice();
//	double var_US_VCS = USMC_VCS.GetVariance();
//
//	USMC_VCS.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_VCS = USMC_VCS.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo and Simples polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_VCS).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo and Simples polynomes  : " << price_US_VCS << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo and Simples polynomes  : " << var_US_VCS << std::endl;
//	std::cout << " US Case - optimal number of simulations for control variable MonteCarlo and Simples polynomes  : " << Optimal_Nb_Simulation_VCS << std::endl;
//
//	delete ugen_VCS;
//	delete ngen_VCS;
//	delete vectorG_VCS;
//	delete BSamerican_VCS;
//
//	////////////////// 2.3.2) Control variable (VDC - Simples Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* vdc_gen_VCS = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_VCS = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_VCS = new NormalBoxMullerVDC(vdc_gen_VCS, vdc_gen_alias_VCS, 0., 1.);
//
//	GaussianVector* vectorG_vdc_VCS = new GaussianVectorCholesky(ngen_vdc_VCS, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_VCS = new BSEulerND(vectorG_vdc_VCS, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_VC_vdc_Simple(N, bsktcall, bsktcallCV, BSamerican_vdc_VCS,
//		basefunc_Simple, price_cf);
//
//	clock::time_point start_US_vdc_VCS = clock::now();
//	USMC_VC_vdc_Simple.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_VCS = clock::now();
//	clock::duration execution_timeUS_vdc_VCS = end_US_vdc_VCS - start_US_vdc_VCS;
//	double price_US_vdc_VCS = USMC_VC_vdc_Simple.GetPrice();
//	double var_US_vdc_VCS = USMC_VC_vdc_Simple.GetVariance();
//
//	USMC_VC_vdc_Simple.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_VCS = USMC_VC_vdc_Simple.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo with VDC sequence and Simple polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_VCS).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo with VDC sequence and Simple polynomes  : " << price_US_vdc_VCS << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo  with VDC sequence and Simple polynomes : " << var_US_vdc_VCS << std::endl;
//	std::cout << " US Case - optimal number of simulations for control variable MonteCarlo with VDC sequence and Simple polynomes : " << Optimal_Nb_Simulation_vdc_VCS << std::endl;
//
//	delete vdc_gen_VCS;
//	delete vdc_gen_alias_VCS;
//	delete ngen_vdc_VCS;
//	delete vectorG_vdc_VCS;
//	delete BSamerican_vdc_VCS;
//
//	////////////////// 2.3.3) Control variable (Sobol -Simple Polynomes n = 2) //////////////////////////////////////
//
//	UniformGenerator* sobol_ugen_CVS = new Sobol();
//	Normal* ngen_sobol_CVS = new NormalBoxMuller(sobol_ugen_CVS, 0., 1.);
//
//	GaussianVector* vectorG_sobol_CVS = new GaussianVectorCholesky(ngen_sobol_CVS, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_CVS = new BSEulerND(vectorG_sobol_CVS, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_sobol_CV_Simple(N, bsktcall, bsktcallCV,BSamerican_sobol_CVS, basefunc_Simple,price_cf);
//
//	clock::time_point start_US_sobol_CVS = clock::now();
//	USMC_sobol_CV_Simple.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_CVS = clock::now();
//	clock::duration execution_timeUS_sobol_CVS = end_US_sobol_CVS - start_US_sobol_CVS;
//	double price_US_sobol_CVS = USMC_sobol_CV_Simple.GetPrice();
//	double var_US_sobol_CVS = USMC_sobol_CV_Simple.GetVariance();
//
//	USMC_sobol_CV_Simple.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_CVS = USMC_sobol_CV_Simple.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo with Sobol sequence and Simple polynomes  : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_CVS).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo with Sobol sequence and Simple polynomes  : " << price_US_sobol_CVS << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo with Sobol sequence and Simple polynomes  : " << var_US_sobol_CVS << std::endl;
//	std::cout << " US Case - optimal number of simulations forcontrol variable  MonteCarlo  with Sobol sequence and Simple polynomes  : " << Optimal_Nb_Simulation_sobol_CVS << std::endl;
//
//	delete sobol_ugen_CVS;
//	delete ngen_sobol_CVS;
//	delete vectorG_sobol_CVS;
//	delete BSamerican_sobol_CVS;
//
//	//////////////////////////// Antithetic /////////////////////////////////////////////////////////////////////
//
//	//////////////////////// 3.1) Antithetic (Ecuyer Combined - Laguerre polynome n = 2) ////////////////////////////
//
//	UniformGenerator* ugen_anti = new EcuyerCombined();
//	Normal* ngen_anti = new NormalBoxMuller(ugen_anti, 0., 1.);
//
//	GaussianVector* vectorG_anti = new GaussianVectorCholesky(ngen_anti, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_anti = new BSEulerNDAntithetic(vectorG_anti, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_anti(N, bsktcall, BSamerican_anti,
//		basefunc_Laguerre);
//
//	clock::time_point start_US_anti = clock::now();
//	USMC_anti.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_anti = clock::now();
//	clock::duration execution_timeUS_anti = end_US_anti - start_US_anti;
//	double price_US_anti = USMC_anti.GetPrice();
//	double var_US_anti = USMC_anti.GetVariance();
//
//	USMC_anti.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_anti = USMC_anti.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo : " << price_US_anti << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo : " << var_US_anti << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo : " << Optimal_Nb_Simulation_anti << std::endl;
//
//	delete ugen_anti;
//	delete ngen_anti;
//	delete vectorG_anti;
//	delete BSamerican_anti;
//
//
//	//////////////////////// 3.1.2) Antithetic (VDC - Laguerre polynome n = 2) ////////////////////////////
//
//	UniformGenerator* vdc_gen_anti = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_anti = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_anti = new NormalBoxMullerVDC(vdc_gen_anti, vdc_gen_alias_anti, 0., 1.);
//
//	GaussianVector* vectorG_vdc_anti = new GaussianVectorCholesky(ngen_vdc_anti, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_anti = new BSEulerNDAntithetic(vectorG_vdc_anti, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_vdc_anti(N, bsktcall, BSamerican_vdc_anti,
//		basefunc_Laguerre);
//
//	clock::time_point start_US_vdc_anti = clock::now();
//	USMC_vdc_anti.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_anti = clock::now();
//	clock::duration execution_timeUS_vdc_anti = end_US_vdc_anti - start_US_vdc_anti;
//	double price_US_vdc_anti = USMC_vdc_anti.GetPrice();
//	double var_US_vdc_anti = USMC_vdc_anti.GetVariance();
//
//	USMC_vdc_anti.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_anti = USMC_vdc_anti.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with VDC sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with VDC sequence : " << price_US_vdc_anti << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with VDC sequence : " << var_US_vdc_anti << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with VDC sequence : " << Optimal_Nb_Simulation_vdc_anti << std::endl;
//
//	delete vdc_gen_anti;
//	delete vdc_gen_alias_anti;
//	delete ngen_vdc_anti;
//	delete vectorG_vdc_anti;
//	delete BSamerican_vdc_anti;
//
//	//////////////////////// 3.1.3) Antithetic (Sobol - Laguerre polynome n = 2) ////////////////////////////
//
//	UniformGenerator* gen_sobol_anti = new Sobol();
//	Normal* ngen_sobol_anti = new NormalBoxMuller(gen_sobol_anti, 0., 1.);
//
//	GaussianVector* vectorG_sobol_anti = new GaussianVectorCholesky(ngen_sobol_anti, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_anti = new BSEulerNDAntithetic(vectorG_sobol_anti, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_sobol_anti(N, bsktcall, BSamerican_sobol_anti,
//		basefunc_Laguerre);
//
//	clock::time_point start_US_sobol_anti = clock::now();
//	USMC_sobol_anti.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_anti = clock::now();
//	clock::duration execution_timeUS_sobol_anti = end_US_sobol_anti - start_US_sobol_anti;
//	double price_US_sobol_anti = USMC_sobol_anti.GetPrice();
//	double var_US_sobol_anti = USMC_sobol_anti.GetVariance();
//
//	USMC_sobol_anti.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_anti = USMC_sobol_anti.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with Sobol sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with Sobol sequence : " << price_US_sobol_anti << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with Sobol sequence : " << var_US_sobol_anti << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with Sobol sequence : " << Optimal_Nb_Simulation_sobol_anti << std::endl;
//
//	delete gen_sobol_anti;
//	delete ngen_sobol_anti;
//	delete vectorG_sobol_anti;
//	delete BSamerican_sobol_anti;
//
//	//////////////////////// 3.2) Antithetic (Ecuyer Combined - Hermite polynome n = 2) ////////////////////////////
//
//	UniformGenerator* ugen_antiH = new EcuyerCombined();
//	Normal* ngen_antiH = new NormalBoxMuller(ugen_antiH, 0., 1.);
//
//	GaussianVector* vectorG_antiH = new GaussianVectorCholesky(ngen_antiH, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_antiH = new BSEulerNDAntithetic(vectorG_antiH, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_anti_Hermite(N, bsktcall, BSamerican_antiH,
//		basefunc_Hermite);
//
//	clock::time_point start_US_antiH = clock::now();
//	USMC_anti_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_antiH = clock::now();
//	clock::duration execution_timeUS_antiH = end_US_antiH - start_US_antiH;
//	double price_US_antiH = USMC_anti_Hermite.GetPrice();
//	double var_US_antiH = USMC_anti_Hermite.GetVariance();
//
//	USMC_anti_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_antiH = USMC_anti_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with Hermite polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_antiH).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with Hermite polynome : " << price_US_antiH << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with Hermite polynome : " << var_US_antiH << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with Hermite polynome : " << Optimal_Nb_Simulation_antiH << std::endl;
//
//	delete ugen_antiH;
//	delete ngen_antiH;
//	delete vectorG_antiH;
//	delete BSamerican_antiH;
//
//
//	//////////////////////// 3.2.2) Antithetic (VDC - Hermite polynome n = 2) ////////////////////////////
//
//	UniformGenerator* vdc_gen_antiH = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_antiH = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_antiH = new NormalBoxMullerVDC(vdc_gen_antiH, vdc_gen_alias_antiH, 0., 1.);
//
//	GaussianVector* vectorG_vdc_antiH = new GaussianVectorCholesky(ngen_vdc_antiH, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_antiH = new BSEulerNDAntithetic(vectorG_vdc_antiH, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_vdc_anti_Hermite(N, bsktcall, BSamerican_vdc_antiH,
//		basefunc_Hermite);
//
//	clock::time_point start_US_vdc_antiH = clock::now();
//	USMC_vdc_anti_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_antiH = clock::now();
//	clock::duration execution_timeUS_vdc_antiH = end_US_vdc_antiH - start_US_vdc_antiH;
//	double price_US_vdc_antiH = USMC_vdc_anti_Hermite.GetPrice();
//	double var_US_vdc_antiH = USMC_vdc_anti_Hermite.GetVariance();
//
//	USMC_vdc_anti_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_antiH = USMC_vdc_anti_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with VDC sequence and Hermite Polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_antiH).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with VDC sequence and Hermite Polynome : " << price_US_vdc_antiH << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with VDC sequence and Hermite Polynome : " << var_US_vdc_antiH << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with VDC sequence and Hermite Polynome : " << Optimal_Nb_Simulation_vdc_antiH << std::endl;
//
//	delete vdc_gen_antiH;
//	delete vdc_gen_alias_antiH;
//	delete ngen_vdc_antiH;
//	delete vectorG_vdc_antiH;
//	delete BSamerican_vdc_antiH;
//
//	//////////////////////// 3.2.3) Antithetic (Sobol - Hermite polynome n = 2) ////////////////////////////
//
//	UniformGenerator* gen_sobol_antiH = new Sobol();
//	Normal* ngen_sobol_antiH = new NormalBoxMuller(gen_sobol_antiH, 0., 1.);
//
//	GaussianVector* vectorG_sobol_antiH = new GaussianVectorCholesky(ngen_sobol_antiH, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_antiH = new BSEulerNDAntithetic(vectorG_sobol_antiH, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_sobol_anti_Hermite(N, bsktcall, BSamerican_sobol_antiH,
//		basefunc_Hermite);
//
//	clock::time_point start_US_sobol_antiH = clock::now();
//	USMC_sobol_anti_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_antiH = clock::now();
//	clock::duration execution_timeUS_sobol_antiH = end_US_sobol_antiH - start_US_sobol_antiH;
//	double price_US_sobol_antiH = USMC_sobol_anti_Hermite.GetPrice();
//	double var_US_sobol_antiH = USMC_sobol_anti_Hermite.GetVariance();
//
//	USMC_sobol_anti_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_antiH = USMC_sobol_anti_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with Sobol sequence and Hermite Polynome: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_antiH).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with Sobol sequence and Hermite Polynome : " << price_US_sobol_antiH << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with Sobol sequence and Hermite Polynome : " << var_US_sobol_antiH << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with Sobol sequence and Hermite Polynome : " << Optimal_Nb_Simulation_sobol_antiH << std::endl;
//
//	delete gen_sobol_antiH;
//	delete ngen_sobol_antiH;
//	delete vectorG_sobol_antiH;
//	delete BSamerican_sobol_antiH;
//
//	//////////////////////// 3.3) Antithetic (Ecuyer Combined - Simple polynome n = 2) ////////////////////////////
//
//	UniformGenerator* ugen_antiS = new EcuyerCombined();
//	Normal* ngen_antiS = new NormalBoxMuller(ugen_antiS, 0., 1.);
//
//	GaussianVector* vectorG_antiS = new GaussianVectorCholesky(ngen_antiS, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_antiS = new BSEulerNDAntithetic(vectorG_antiS, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_anti_Simple(N, bsktcall, BSamerican_antiS,
//		basefunc_Simple);
//
//	clock::time_point start_US_antiS = clock::now();
//	USMC_anti_Simple.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_antiS = clock::now();
//	clock::duration execution_timeUS_antiS = end_US_antiS - start_US_antiS;
//	double price_US_antiS = USMC_anti_Simple.GetPrice();
//	double var_US_antiS = USMC_anti_Simple.GetVariance();
//
//	USMC_anti_Simple.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_antiS = USMC_anti_Simple.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with Simple polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_antiS).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with Simple polynome : " << price_US_antiS << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with Simple polynome : " << var_US_antiS << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with Simple polynome : " << Optimal_Nb_Simulation_antiS << std::endl;
//
//	delete ugen_antiS;
//	delete ngen_antiS;
//	delete vectorG_antiS;
//	delete BSamerican_antiS;
//
//
//	//////////////////////// 3.3.2) Antithetic (VDC - Simple polynome n = 2) ////////////////////////////
//
//	UniformGenerator* vdc_gen_antiS = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_antiS = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_antiS = new NormalBoxMullerVDC(vdc_gen_antiS, vdc_gen_alias_antiS, 0., 1.);
//
//	GaussianVector* vectorG_vdc_antiS = new GaussianVectorCholesky(ngen_vdc_antiS, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_antiS = new BSEulerNDAntithetic(vectorG_vdc_antiS, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_vdc_anti_Simple(N, bsktcall, BSamerican_vdc_antiS,
//		basefunc_Simple);
//
//	clock::time_point start_US_vdc_antiS = clock::now();
//	USMC_vdc_anti_Simple.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_antiS = clock::now();
//	clock::duration execution_timeUS_vdc_antiS = end_US_vdc_antiS - start_US_vdc_antiS;
//	double price_US_vdc_antiS = USMC_vdc_anti_Simple.GetPrice();
//	double var_US_vdc_antiS = USMC_vdc_anti_Simple.GetVariance();
//
//	USMC_vdc_anti_Simple.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_antiS = USMC_vdc_anti_Simple.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with VDC sequence and Simple Polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_antiS).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with VDC sequence and Simple Polynome : " << price_US_vdc_antiS << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with VDC sequence and Simple Polynome : " << var_US_vdc_antiS << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with VDC sequence and Simple Polynome : " << Optimal_Nb_Simulation_vdc_antiS << std::endl;
//
//	delete vdc_gen_antiS;
//	delete vdc_gen_alias_antiS;
//	delete ngen_vdc_antiS;
//	delete vectorG_vdc_antiS;
//	delete BSamerican_vdc_antiS;
//
//	//////////////////////// 3.3.3) Antithetic (Sobol - Simple polynome n = 2) ////////////////////////////
//
//	UniformGenerator* gen_sobol_antiS = new Sobol();
//	Normal* ngen_sobol_antiS = new NormalBoxMuller(gen_sobol_antiS, 0., 1.);
//
//	GaussianVector* vectorG_sobol_antiS = new GaussianVectorCholesky(ngen_sobol_antiS, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_antiS = new BSEulerNDAntithetic(vectorG_sobol_antiS, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_sobol_anti_Simple(N, bsktcall, BSamerican_sobol_antiS,
//		basefunc_Simple);
//
//	clock::time_point start_US_sobol_antiS = clock::now();
//	USMC_sobol_anti_Simple.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_antiS = clock::now();
//	clock::duration execution_timeUS_sobol_antiS = end_US_sobol_antiS - start_US_sobol_antiS;
//	double price_US_sobol_antiS = USMC_sobol_anti_Simple.GetPrice();
//	double var_US_sobol_antiS = USMC_sobol_anti_Simple.GetVariance();
//
//	USMC_sobol_anti_Simple.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_antiS = USMC_sobol_anti_Simple.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with Sobol sequence and Simple Polynome: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_antiS).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with Sobol sequence and Simple Polynome : " << price_US_sobol_antiS << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with Sobol sequence and Simple Polynome : " << var_US_sobol_antiS << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with Sobol sequence and Simple Polynome : " << Optimal_Nb_Simulation_sobol_antiS << std::endl;
//
//	delete gen_sobol_antiS;
//	delete ngen_sobol_antiS;
//	delete vectorG_sobol_antiS;
//	delete BSamerican_sobol_antiS;
//
//	//////////////////////////// VC + Antithetic /////////////////////////////////////////////////////////////////////
//
////////////////////////// 4.1) VC + Antithetic (Ecuyer Combined - Laguerre polynome n = 2) ////////////////////////////
//
//	UniformGenerator* ugen_anti_cv = new EcuyerCombined();
//	Normal* ngen_anti_cv = new NormalBoxMuller(ugen_anti_cv, 0., 1.);
//
//	GaussianVector* vectorG_anti_cv = new GaussianVectorCholesky(ngen_anti_cv, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_anti_cv = new BSEulerNDAntithetic(vectorG_anti_cv, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_anti_CV(N, bsktcall, bsktcallCV, BSamerican_anti_cv, basefunc_Laguerre, price_cf);
//
//
//	clock::time_point start_US_anti_cv = clock::now();
//	USMC_anti_CV.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_anti_cv = clock::now();
//	clock::duration execution_timeUS_anti_cv = end_US_anti_cv - start_US_anti_cv;
//	double price_US_anti_cv = USMC_anti_CV.GetPrice();
//	double var_US_anti_cv = USMC_anti_CV.GetVariance();
//
//	USMC_anti_CV.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_anti_cv = USMC_anti_CV.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti_cv).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo : " << price_US_anti_cv << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo : " << var_US_anti_cv << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo : " << Optimal_Nb_Simulation_anti_cv << std::endl;
//
//	delete ugen_anti_cv;
//	delete ngen_anti_cv;
//	delete vectorG_anti_cv;
//	delete BSamerican_anti_cv;
//
//	//////////////////////// 4.1.2) VC + Antithetic (VDC - Laguerre polynome n = 2) ////////////////////////////
//
//	UniformGenerator* vdc_gen_anti_cv = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_anti_cv = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_anti_cv = new NormalBoxMullerVDC(vdc_gen_anti_cv, vdc_gen_alias_anti_cv, 0., 1.);
//
//	GaussianVector* vectorG_vdc_anti_cv = new GaussianVectorCholesky(ngen_vdc_anti_cv, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_anti_cv = new BSEulerNDAntithetic(vectorG_vdc_anti_cv, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_vdc_anti_CV(N, bsktcall, bsktcallCV, BSamerican_vdc_anti_cv, basefunc_Laguerre, price_cf);
//
//
//	clock::time_point start_US_vdc_anti_cv = clock::now();
//	USMC_vdc_anti_CV.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_anti_cv = clock::now();
//	clock::duration execution_timeUS_vdc_anti_cv = end_US_vdc_anti_cv - start_US_vdc_anti_cv;
//	double price_US_vdc_anti_cv = USMC_vdc_anti_CV.GetPrice();
//	double var_US_vdc_anti_cv = USMC_vdc_anti_CV.GetVariance();
//
//	USMC_vdc_anti_CV.OptimalNbSimul(tolerated_error);
//	size_t Optimal_vdc_Nb_Simulation_anti_cv = USMC_vdc_anti_CV.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo with VDC sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti_cv).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo  with VDC sequence: " << price_US_vdc_anti_cv << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo with VDC sequence : " << var_US_vdc_anti_cv << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo with VDC sequence : " << Optimal_vdc_Nb_Simulation_anti_cv << std::endl;
//
//	delete vdc_gen_anti_cv;
//	delete vdc_gen_alias_anti_cv;
//	delete ngen_vdc_anti_cv;
//	delete vectorG_vdc_anti_cv;
//	delete BSamerican_vdc_anti_cv;
//
//	//////////////////////// 4.1.2) VC + Antithetic (Sobol - Laguerre polynome n = 2) ////////////////////////////
//
//	UniformGenerator* sobol_gen_anti_cv = new Sobol();
//	Normal* ngen_sobol_anti_cv = new NormalBoxMuller(sobol_gen_anti_cv, 0., 1.);
//
//	GaussianVector* vectorG_sobol_anti_cv = new GaussianVectorCholesky(ngen_sobol_anti_cv, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_anti_cv = new BSEulerNDAntithetic(vectorG_sobol_anti_cv, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_sobol_anti_CV(N, bsktcall, bsktcallCV, BSamerican_sobol_anti_cv, basefunc_Laguerre, price_cf);
//
//
//	clock::time_point start_US_sobol_anti_cv = clock::now();
//	USMC_sobol_anti_CV.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_anti_cv = clock::now();
//	clock::duration execution_timeUS_sobol_anti_cv = end_US_sobol_anti_cv - start_US_sobol_anti_cv;
//	double price_US_sobol_anti_cv = USMC_sobol_anti_CV.GetPrice();
//	double var_US_sobol_anti_cv = USMC_sobol_anti_CV.GetVariance();
//
//	USMC_sobol_anti_CV.OptimalNbSimul(tolerated_error);
//	size_t Optimal_sobol_Nb_Simulation_anti_cv = USMC_sobol_anti_CV.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo with Sobol sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti_cv).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo  with Sobol sequence: " << price_US_sobol_anti_cv << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo with Sobol sequence : " << var_US_sobol_anti_cv << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo with Sobol sequence : " << Optimal_sobol_Nb_Simulation_anti_cv << std::endl;
//
//	delete sobol_gen_anti_cv;
//	delete ngen_sobol_anti_cv;
//	delete vectorG_sobol_anti_cv;
//	delete BSamerican_sobol_anti_cv;
//
//	//////////////////////// 4.2) VC + Antithetic (Ecuyer Combined - Hermite polynome n = 2) ////////////////////////////
//
//	UniformGenerator* ugen_anti_cvH = new EcuyerCombined();
//	Normal* ngen_anti_cvH = new NormalBoxMuller(ugen_anti_cvH, 0., 1.);
//
//	GaussianVector* vectorG_anti_cvH = new GaussianVectorCholesky(ngen_anti_cvH, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_anti_cvH = new BSEulerNDAntithetic(vectorG_anti_cvH, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_anti_CV_Hermite(N, bsktcall, bsktcallCV, BSamerican_anti_cvH, basefunc_Hermite, price_cf);
//
//
//	clock::time_point start_US_anti_cvH = clock::now();
//	USMC_anti_CV_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_anti_cvH = clock::now();
//	clock::duration execution_timeUS_anti_cvH = end_US_anti_cvH - start_US_anti_cvH;
//	double price_US_anti_cvH = USMC_anti_CV_Hermite.GetPrice();
//	double var_US_anti_cvH = USMC_anti_CV_Hermite.GetVariance();
//
//	USMC_anti_CV_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_anti_cvH = USMC_anti_CV_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo and Hermite polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti_cvH).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo and Hermite polynome : " << price_US_anti_cvH << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo and Hermite polynome : " << var_US_anti_cvH << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo and Hermite Polynome : " << Optimal_Nb_Simulation_anti_cvH << std::endl;
//
//	delete ugen_anti_cvH;
//	delete ngen_anti_cvH;
//	delete vectorG_anti_cvH;
//	delete BSamerican_anti_cvH;
//
//	//////////////////////// 4.2.2) VC + Antithetic (VDC - Hermite polynome n = 2) ////////////////////////////
//
//	UniformGenerator* vdc_gen_anti_cvH = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_anti_cvH = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_anti_cvH = new NormalBoxMullerVDC(vdc_gen_anti_cvH, vdc_gen_alias_anti_cvH, 0., 1.);
//
//	GaussianVector* vectorG_vdc_anti_cvH = new GaussianVectorCholesky(ngen_vdc_anti_cvH, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_anti_cvH = new BSEulerNDAntithetic(vectorG_vdc_anti_cvH, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_vdc_anti_CV_Hermite(N, bsktcall, bsktcallCV, BSamerican_vdc_anti_cvH, basefunc_Hermite, price_cf);
//
//
//	clock::time_point start_US_vdc_anti_cvH = clock::now();
//	USMC_vdc_anti_CV_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_anti_cvH = clock::now();
//	clock::duration execution_timeUS_vdc_anti_cvH = end_US_vdc_anti_cvH - start_US_vdc_anti_cvH;
//	double price_US_vdc_anti_cvH = USMC_vdc_anti_CV_Hermite.GetPrice();
//	double var_US_vdc_anti_cvH = USMC_vdc_anti_CV_Hermite.GetVariance();
//
//	USMC_vdc_anti_CV_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_vdc_Nb_Simulation_anti_cvH = USMC_vdc_anti_CV_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo with VDC sequence and Hermite polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti_cvH).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo  with VDC sequence and Hermite polynome : " << price_US_vdc_anti_cvH << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo with VDC sequence and Hermite polynome : " << var_US_vdc_anti_cvH << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo with VDC sequence and Hermite polynome : " << Optimal_vdc_Nb_Simulation_anti_cvH << std::endl;
//
//	delete vdc_gen_anti_cvH;
//	delete vdc_gen_alias_anti_cvH;
//	delete ngen_vdc_anti_cvH;
//	delete vectorG_vdc_anti_cvH;
//	delete BSamerican_vdc_anti_cvH;
//
//	//////////////////////// 4.2.3) VC + Antithetic (Sobol - Hermite polynome n = 2) ////////////////////////////
//
//	UniformGenerator* sobol_gen_anti_cvH = new Sobol();
//	Normal* ngen_sobol_anti_cvH = new NormalBoxMuller(sobol_gen_anti_cvH, 0., 1.);
//
//	GaussianVector* vectorG_sobol_anti_cvH = new GaussianVectorCholesky(ngen_sobol_anti_cvH, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_anti_cvH = new BSEulerNDAntithetic(vectorG_sobol_anti_cvH, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_sobol_anti_CV_Hermite(N, bsktcall, bsktcallCV, BSamerican_sobol_anti_cvH, basefunc_Hermite, price_cf);
//
//
//	clock::time_point start_US_sobol_anti_cvH = clock::now();
//	USMC_sobol_anti_CV_Hermite.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_anti_cvH = clock::now();
//	clock::duration execution_timeUS_sobol_anti_cvH = end_US_sobol_anti_cvH - start_US_sobol_anti_cvH;
//	double price_US_sobol_anti_cvH = USMC_sobol_anti_CV_Hermite.GetPrice();
//	double var_US_sobol_anti_cvH = USMC_sobol_anti_CV_Hermite.GetVariance();
//
//	USMC_sobol_anti_CV_Hermite.OptimalNbSimul(tolerated_error);
//	size_t Optimal_sobol_Nb_Simulation_anti_cvH = USMC_sobol_anti_CV_Hermite.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo with Sobol sequence and Hermite polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti_cvH).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo  with Sobol sequence and Hermite polynome : " << price_US_sobol_anti_cvH << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo with Sobol sequence and Hermite polynome : " << var_US_sobol_anti_cvH << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo with Sobol sequence and Hermite polynome : " << Optimal_sobol_Nb_Simulation_anti_cvH << std::endl;
//
//	delete sobol_gen_anti_cvH;
//	delete ngen_sobol_anti_cvH;
//	delete vectorG_sobol_anti_cvH;
//	delete BSamerican_sobol_anti_cvH;
//
//	//////////////////////// 4.3) VC + Antithetic (Ecuyer Combined - simple polynome n = 2) ////////////////////////////
//
//	UniformGenerator* ugen_anti_cvS = new EcuyerCombined();
//	Normal* ngen_anti_cvS = new NormalBoxMuller(ugen_anti_cvS, 0., 1.);
//
//	GaussianVector* vectorG_anti_cvS = new GaussianVectorCholesky(ngen_anti_cvS, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_anti_cvS = new BSEulerNDAntithetic(vectorG_anti_cvS, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_anti_CV_Simple(N, bsktcall, bsktcallCV, BSamerican_anti_cvS, basefunc_Simple, price_cf);
//
//
//	clock::time_point start_US_anti_cvS = clock::now();
//	USMC_anti_CV_Simple.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_anti_cvS = clock::now();
//	clock::duration execution_timeUS_anti_cvS = end_US_anti_cvS - start_US_anti_cvS;
//	double price_US_anti_cvS = USMC_anti_CV_Simple.GetPrice();
//	double var_US_anti_cvS = USMC_anti_CV_Simple.GetVariance();
//
//	USMC_anti_CV_Simple.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_anti_cvS = USMC_anti_CV_Simple.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo and Simple polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti_cvS).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo and Simple polynome : " << price_US_anti_cvS << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo and Simple polynome : " << var_US_anti_cvS << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo and Simple Polynome : " << Optimal_Nb_Simulation_anti_cvS << std::endl;
//
//	delete ugen_anti_cvS;
//	delete ngen_anti_cvS;
//	delete vectorG_anti_cvS;
//	delete BSamerican_anti_cvS;
//
//	//////////////////////// 4.3.2) VC + Antithetic (VDC - Simple polynome n = 2) ////////////////////////////
//
//	UniformGenerator* vdc_gen_anti_cvS = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_anti_cvS = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_anti_cvS = new NormalBoxMullerVDC(vdc_gen_anti_cvS, vdc_gen_alias_anti_cvS, 0., 1.);
//
//	GaussianVector* vectorG_vdc_anti_cvS = new GaussianVectorCholesky(ngen_vdc_anti_cvS, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_anti_cvS = new BSEulerNDAntithetic(vectorG_vdc_anti_cvS, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_vdc_anti_CV_Simple(N, bsktcall, bsktcallCV, BSamerican_vdc_anti_cvS, basefunc_Simple, price_cf);
//
//
//	clock::time_point start_US_vdc_anti_cvS = clock::now();
//	USMC_vdc_anti_CV_Simple.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_anti_cvS = clock::now();
//	clock::duration execution_timeUS_vdc_anti_cvS = end_US_vdc_anti_cvS - start_US_vdc_anti_cvS;
//	double price_US_vdc_anti_cvS = USMC_vdc_anti_CV_Simple.GetPrice();
//	double var_US_vdc_anti_cvS = USMC_vdc_anti_CV_Simple.GetVariance();
//
//	USMC_vdc_anti_CV_Simple.OptimalNbSimul(tolerated_error);
//	size_t Optimal_vdc_Nb_Simulation_anti_cvS = USMC_vdc_anti_CV_Simple.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo with VDC sequence and simple polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti_cvS).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo  with VDC sequence and simple polynome : " << price_US_vdc_anti_cvS << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo with VDC sequence and simple polynome : " << var_US_vdc_anti_cvS << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo with VDC sequence and simple polynome : " << Optimal_vdc_Nb_Simulation_anti_cvS << std::endl;
//
//	delete vdc_gen_anti_cvS;
//	delete vdc_gen_alias_anti_cvS;
//	delete ngen_vdc_anti_cvS;
//	delete vectorG_vdc_anti_cvS;
//	delete BSamerican_vdc_anti_cvS;
//
//	//////////////////////// 4.3.3) VC + Antithetic (Sobol - simple polynome n = 2) ////////////////////////////
//
//	UniformGenerator* sobol_gen_anti_cvS = new Sobol();
//	Normal* ngen_sobol_anti_cvS = new NormalBoxMuller(sobol_gen_anti_cvS, 0., 1.);
//
//	GaussianVector* vectorG_sobol_anti_cvS = new GaussianVectorCholesky(ngen_sobol_anti_cvS, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_anti_cvS = new BSEulerNDAntithetic(vectorG_sobol_anti_cvS, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_sobol_anti_CV_Simple(N, bsktcall, bsktcallCV, BSamerican_sobol_anti_cvS, basefunc_Simple, price_cf);
//
//
//	clock::time_point start_US_sobol_anti_cvS = clock::now();
//	USMC_sobol_anti_CV_Simple.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_anti_cvS = clock::now();
//	clock::duration execution_timeUS_sobol_anti_cvS = end_US_sobol_anti_cvS - start_US_sobol_anti_cvS;
//	double price_US_sobol_anti_cvS = USMC_sobol_anti_CV_Simple.GetPrice();
//	double var_US_sobol_anti_cvS = USMC_sobol_anti_CV_Simple.GetVariance();
//
//	USMC_sobol_anti_CV_Simple.OptimalNbSimul(tolerated_error);
//	size_t Optimal_sobol_Nb_Simulation_anti_cvS = USMC_sobol_anti_CV_Simple.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo with Sobol sequence and simple polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti_cvS).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo  with Sobol sequence and simple polynome : " << price_US_sobol_anti_cvS << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo with Sobol sequence and Simple polynome : " << var_US_sobol_anti_cvS << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo with Sobol sequence and Simple polynome : " << Optimal_sobol_Nb_Simulation_anti_cvS << std::endl;
//
//	delete sobol_gen_anti_cvS;
//	delete ngen_sobol_anti_cvS;
//	delete vectorG_sobol_anti_cvS;
//	delete BSamerican_sobol_anti_cvS;

///// BERMUDAN PART //////////////////////////////////////////////////////////////////////////////

	//creattion of the exercice schedule for the Bermudan option

size_t n = 12;
matrix Schedule_exec(n, 1);

double test = 0.08;

for (size_t i = 0; i < n; i++)
{

	Schedule_exec(i, 0) = test;
	test = 0.08 * (i + 1.);
}


//creation of the calendar management to get th right spot vector at each simulation 

CalendarManagement* wkday = new rounded_workingdays(0);

/////////////////// NO VARIANCE REDUCTION /////////////////////////////

////////////////// 1.1) Vanille (Ecuyer Combined - Laguerre Polynomes n = 2) //////////////////////////////////////

UniformGenerator* ugen_berm = new EcuyerCombined();
Normal* ngen_berm = new NormalBoxMuller(ugen, 0., 1.);

GaussianVector* vectorGberm = new GaussianVectorCholesky(ngen_berm, Sigma, Correl, CovarMatrix);
RandomProcess* BSbermudan = new BSEulerND(vectorGberm, spot_m, rate);

Bermudean_BasketOption berm(N, bsktcall, BSbermudan, basefunc_Laguerre, wkday, Schedule_exec);

clock::time_point start_berm = clock::now();
berm.Simulate(startTime, endTime, nbsteps);
clock::time_point end_berm = clock::now();
clock::duration execution_timeberm = end_berm - start_berm;
double price_berm = berm.GetPrice();
double var_berm = berm.GetVariance();

berm.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation = berm.GetNbSimul();

std::cout << " berm Case - exec time for vanilla MonteCarlo : " << std::chrono::duration <double, std::ratio<1>>(execution_timeberm).count() << std::endl;
std::cout << " berm Case - price for vanilla MonteCarlo : " << price_berm << std::endl;
std::cout << " berm Case - variance for vanilla MonteCarlo : " << var_berm << std::endl;
std::cout << " berm Case - optimal number of simulations for vanilla MonteCarlo : " << Optimal_Nb_Simulation << std::endl;

delete ugen_berm;
delete ngen_berm;
delete vectorGberm;
delete BSbermudan;

/*
////////////////// 1.1.2) Vanille (VDC - Laguerre Polynomes n = 2) //////////////////////////////////////

UniformGenerator* vdc_gen = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias = new VanDerCorput(11, 1);
Normal* ngen_vdc = new NormalBoxMullerVDC(vdc_gen, vdc_gen_alias, 0., 1.);

GaussianVector* vectorG_vdc = new GaussianVectorCholesky(ngen_vdc, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc = new BSEulerND(vectorG_vdc, spot_m, rate);

Bermudean_BasketOption USMC_vdc(N, bsktcall, BSamerican_vdc, basefunc_Laguerre, wkday, Schedule_exec);

clock::time_point start_US_vdc = clock::now();
USMC_vdc.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc = clock::now();
clock::duration execution_timeUS_vdc = end_US_vdc - start_US_vdc;
double price_US_vdc = USMC_vdc.GetPrice();
double var_US_vdc = USMC_vdc.GetVariance();

USMC_vdc.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_vdc = USMC_vdc.GetNbSimul();

std::cout << " Berm Case - exec time for vanilla MonteCarlo with VDC sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc).count() << std::endl;
std::cout << " Berm Case - price for vanilla MonteCarlo with VDC sequence : " << price_US_vdc << std::endl;
std::cout << " Berm Case - variance for vanilla MonteCarlo  with VDC sequence : " << var_US_vdc << std::endl;
std::cout << " Berm Case - optimal number of simulations for vanilla MonteCarlo  with VDC sequence : " << Optimal_Nb_Simulation_vdc << std::endl;

delete vdc_gen;
delete vdc_gen_alias;
delete ngen_vdc;
delete vectorG_vdc;
delete BSamerican_vdc;

////////////////// 1.1.3) Vanille (Sobol - Laguerre Polynomes n = 2) //////////////////////////////////////

UniformGenerator* sobol_ugen = new Sobol();
Normal* ngen_sobol = new NormalBoxMuller(sobol_ugen, 0., 1.);

GaussianVector* vectorG_sobol = new GaussianVectorCholesky(ngen_sobol, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol = new BSEulerND(vectorG_sobol, spot_m, rate);

Bermudean_BasketOption USMC_sobol(N, bsktcall, BSamerican_sobol, basefunc_Laguerre, wkday, Schedule_exec);

clock::time_point start_US_sobol = clock::now();
USMC_sobol.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol = clock::now();
clock::duration execution_timeUS_sobol = end_US_sobol - start_US_sobol;
double price_US_sobol = USMC_sobol.GetPrice();
double var_US_sobol = USMC_sobol.GetVariance();

USMC_sobol.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_sobol = USMC_sobol.GetNbSimul();

std::cout << " Berm Case - exec time for vanilla MonteCarlo with Sobol sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol).count() << std::endl;
std::cout << " Berm Case - price for vanilla MonteCarlo with Sobol sequence : " << price_US_sobol << std::endl;
std::cout << " Berm Case - variance for vanilla MonteCarlo  with Sobol sequence : " << var_US_sobol << std::endl;
std::cout << " Berm Case - optimal number of simulations for vanilla MonteCarlo  with Sobol sequence : " << Optimal_Nb_Simulation_sobol << std::endl;

delete sobol_ugen;
delete ngen_sobol;
delete vectorG_sobol;
delete BSamerican_sobol;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// 1.2) Vanille (Ecuyer Combined - Hermite Polynomes n = 2) //////////////////////////////////////

UniformGenerator* ugen_H = new EcuyerCombined();
Normal* ngen_H = new NormalBoxMuller(ugen_H, 0., 1.);

GaussianVector* vectorG_H = new GaussianVectorCholesky(ngen_H, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_H = new BSEulerND(vectorG_H, spot_m, rate);

Bermudean_BasketOption USMC_Hermite(N, bsktcall, BSamerican_H, basefunc_Hermite, wkday, Schedule_exec);

clock::time_point start_US_H = clock::now();
USMC_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_H = clock::now();
clock::duration execution_timeUS_H = end_US_H - start_US_H;
double price_US_H = USMC_Hermite.GetPrice();
double var_US_H = USMC_Hermite.GetVariance();

USMC_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_H = USMC_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for vanilla MonteCarlo with Hermite polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_H).count() << std::endl;
std::cout << " Berm Case - price for vanilla MonteCarlo with Hermite polynome : " << price_US_H << std::endl;
std::cout << " Berm Case - variance for vanilla MonteCarlo with Hermite polynome : " << var_US_H << std::endl;
std::cout << " Berm Case - optimal number of simulations for vanilla MonteCarlo with Hermite polynome : " << Optimal_Nb_Simulation_H << std::endl;

delete ugen_H;
delete ngen_H;
delete vectorG_H;
delete BSamerican_H;

////////////////// 1.2.2) Vanille (VDC - Hermite Polynomes n = 2) //////////////////////////////////////

UniformGenerator* vdc_gen_H = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_H = new VanDerCorput(11, 1);
Normal* ngen_vdc_H = new NormalBoxMullerVDC(vdc_gen_H, vdc_gen_alias_H, 0., 1.);

GaussianVector* vectorG_vdc_H = new GaussianVectorCholesky(ngen_vdc_H, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_H = new BSEulerND(vectorG_vdc_H, spot_m, rate);

Bermudean_BasketOption USMC_vdc_Hermite(N, bsktcall, BSamerican_vdc_H, basefunc_Hermite, wkday, Schedule_exec);

clock::time_point start_US_vdc_H = clock::now();
USMC_vdc_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_H = clock::now();
clock::duration execution_timeUS_vdc_H = end_US_vdc_H - start_US_vdc_H;
double price_US_vdc_H = USMC_vdc_Hermite.GetPrice();
double var_US_vdc_H = USMC_vdc_Hermite.GetVariance();

USMC_vdc_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_vdc_Hermite = USMC_vdc_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for vanilla MonteCarlo with VDC sequence and Hermite Polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_H).count() << std::endl;
std::cout << " Berm Case - price for vanilla MonteCarlo with VDC sequence and Hermite Polynomes : " << price_US_vdc_H << std::endl;
std::cout << " Berm Case - variance for vanilla MonteCarlo  with VDC sequence and Hermite Polynomes : " << var_US_vdc_H << std::endl;
std::cout << " Berm Case - optimal number of simulations for vanilla MonteCarlo  with VDC sequence and Hermite Polynomes : " << Optimal_Nb_Simulation_vdc_Hermite << std::endl;

delete vdc_gen_H;
delete vdc_gen_alias_H;
delete ngen_vdc_H;
delete vectorG_vdc_H;
delete BSamerican_vdc_H;

////////////////// 1.2.3) Vanille (Sobol - Hermite Polynomes n = 2) //////////////////////////////////////

UniformGenerator* sobol_ugen_H = new Sobol();
Normal* ngen_sobol_H = new NormalBoxMuller(sobol_ugen_H, 0., 1.);

GaussianVector* vectorG_sobol_H = new GaussianVectorCholesky(ngen_sobol_H, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_H = new BSEulerND(vectorG_sobol_H, spot_m, rate);

Bermudean_BasketOption USMC_sobol_Hermite(N, bsktcall, BSamerican_sobol_H, basefunc_Hermite, wkday, Schedule_exec);

clock::time_point start_US_sobol_H = clock::now();
USMC_sobol_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_H = clock::now();
clock::duration execution_timeUS_sobol_H = end_US_sobol_H - start_US_sobol_H;
double price_US_sobol_H = USMC_sobol_Hermite.GetPrice();
double var_US_sobol_H = USMC_sobol_Hermite.GetVariance();

USMC_sobol_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_sobol_H = USMC_sobol_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for vanilla MonteCarlo with Sobol sequence and Hermite Polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_H).count() << std::endl;
std::cout << " Berm Case - price for vanilla MonteCarlo with Sobol sequence and Hermite Polynomes : " << price_US_sobol_H << std::endl;
std::cout << " Berm Case - variance for vanilla MonteCarlo  with Sobol sequence and Hermite Polynomes : " << var_US_sobol_H << std::endl;
std::cout << " Berm Case - optimal number of simulations for vanilla MonteCarlo  with Sobol sequence and Hermite Polynomes : " << Optimal_Nb_Simulation_sobol_H << std::endl;

delete sobol_ugen_H;
delete ngen_sobol_H;
delete vectorG_sobol_H;
delete BSamerican_sobol_H;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// 1.3) Vanille (Ecuyer Combined - Simples Polynomes n = 2) //////////////////////////////////////

UniformGenerator* ugen_S = new EcuyerCombined();
Normal* ngen_S = new NormalBoxMuller(ugen_S, 0., 1.);

GaussianVector* vectorG_S = new GaussianVectorCholesky(ngen_S, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_S = new BSEulerND(vectorG_S, spot_m, rate);

Bermudean_BasketOption USMC_Simples(N, bsktcall, BSamerican_S, basefunc_Simple, wkday, Schedule_exec);

clock::time_point start_US_S = clock::now();
USMC_Simples.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_S = clock::now();
clock::duration execution_timeUS_S = end_US_S - start_US_S;
double price_US_S = USMC_Simples.GetPrice();
double var_US_S = USMC_Simples.GetVariance();

USMC_Simples.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_S = USMC_Simples.GetNbSimul();

std::cout << " Berm Case - exec time for vanilla MonteCarlo with simples polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_S).count() << std::endl;
std::cout << " Berm Case - price for vanilla MonteCarlo with simples polynome : " << price_US_S << std::endl;
std::cout << " Berm Case - variance for vanilla MonteCarlo with simples polynome : " << var_US_S << std::endl;
std::cout << " Berm Case - optimal number of simulations for vanilla MonteCarlo with simples polynome : " << Optimal_Nb_Simulation_S << std::endl;

delete ugen_S;
delete ngen_S;
delete vectorG_S;
delete BSamerican_S;

////////////////// 1.3.2) Vanille (VDC - simples Polynomes n = 2) //////////////////////////////////////

UniformGenerator* vdc_gen_S = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_S = new VanDerCorput(11, 1);
Normal* ngen_vdc_S = new NormalBoxMullerVDC(vdc_gen_S, vdc_gen_alias_S, 0., 1.);

GaussianVector* vectorG_vdc_S = new GaussianVectorCholesky(ngen_vdc_S, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_S = new BSEulerND(vectorG_vdc_S, spot_m, rate);

Bermudean_BasketOption USMC_vdc_Simples(N, bsktcall, BSamerican_vdc_S, basefunc_Simple, wkday, Schedule_exec);

clock::time_point start_US_vdc_S = clock::now();
USMC_vdc_Simples.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_S = clock::now();
clock::duration execution_timeUS_vdc_S = end_US_vdc_S - start_US_vdc_S;
double price_US_vdc_S = USMC_vdc_Simples.GetPrice();
double var_US_vdc_S = USMC_vdc_Simples.GetVariance();

USMC_vdc_Simples.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_vdc_Simples = USMC_vdc_Simples.GetNbSimul();

std::cout << " Berm Case - exec time for vanilla MonteCarlo with VDC sequence and simples Polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_S).count() << std::endl;
std::cout << " Berm Case - price for vanilla MonteCarlo with VDC sequence and simples Polynomes : " << price_US_vdc_S << std::endl;
std::cout << " Berm Case - variance for vanilla MonteCarlo  with VDC sequence and simples Polynomes : " << var_US_vdc_S << std::endl;
std::cout << " Berm Case - optimal number of simulations for vanilla MonteCarlo  with VDC sequence and simples Polynomes : " << Optimal_Nb_Simulation_vdc_Simples << std::endl;

delete vdc_gen_S;
delete vdc_gen_alias_S;
delete ngen_vdc_S;
delete vectorG_vdc_S;
delete BSamerican_vdc_S;

////////////////// 1.3.3) Vanille (Sobol - simples Polynomes n = 2) //////////////////////////////////////

UniformGenerator* sobol_ugen_S = new Sobol();
Normal* ngen_sobol_S = new NormalBoxMuller(sobol_ugen_S, 0., 1.);

GaussianVector* vectorG_sobol_S = new GaussianVectorCholesky(ngen_sobol_S, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_S = new BSEulerND(vectorG_sobol_S, spot_m, rate);

Bermudean_BasketOption USMC_sobol_Simples(N, bsktcall, BSamerican_sobol_S, basefunc_Simple, wkday, Schedule_exec);

clock::time_point start_US_sobol_S = clock::now();
USMC_sobol_Simples.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_S = clock::now();
clock::duration execution_timeUS_sobol_S = end_US_sobol_S - start_US_sobol_S;
double price_US_sobol_S = USMC_sobol_Simples.GetPrice();
double var_US_sobol_S = USMC_sobol_Simples.GetVariance();

USMC_sobol_Simples.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_sobol_S = USMC_sobol_Simples.GetNbSimul();

std::cout << " Berm Case - exec time for vanilla MonteCarlo with Sobol sequence and simples Polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_S).count() << std::endl;
std::cout << " Berm Case - price for vanilla MonteCarlo with Sobol sequence and simples Polynomes : " << price_US_sobol_S << std::endl;
std::cout << " Berm Case - variance for vanilla MonteCarlo  with Sobol sequence and simples Polynomes : " << var_US_sobol_S << std::endl;
std::cout << " Berm Case - optimal number of simulations for vanilla MonteCarlo  with Sobol sequence and simples Polynomes : " << Optimal_Nb_Simulation_sobol_S << std::endl;

delete sobol_ugen_S;
delete ngen_sobol_S;
delete vectorG_sobol_S;
delete BSamerican_sobol_S;

///////////////////////////////////// VARIANCE REDUCTION TECHNIQUES ////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////// CONTROL VARIABLE /////////////////////////////////////////////////////////////////////


	///create both the control variate payoff and its close_formula price in order to proceed to the contrl variate technique

PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(W, spot_m, 100.);

ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(W, spot_m, CovarMatrix, K,
	rate, endTime);

double df = exp(-rate * endTime);
double price_cf = CFbaskt->operator()(spot_m, df);

//////////////////////// 2.1) Control Variable (Ecuyer Combined - Laguerre polynome n = 2) ////////////////////////////


UniformGenerator* ugen_VC = new EcuyerCombined();
Normal* ngen_VC = new NormalBoxMuller(ugen_VC, 0., 1.);

GaussianVector* vectorG_VC = new GaussianVectorCholesky(ngen_VC, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_VC = new BSEulerND(vectorG_VC, spot_m, rate);

Bermudean_BasketOption_CV USMC_VC(N, bsktcall, bsktcallCV, BSamerican_VC,
	basefunc_Laguerre, wkday, Schedule_exec, price_cf);

clock::time_point start_US_VC = clock::now();
USMC_VC.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_VC = clock::now();
clock::duration execution_timeUS_VC = end_US_VC - start_US_VC;
double price_US_VC = USMC_VC.GetPrice();
double var_US_VC = USMC_VC.GetVariance();

USMC_VC.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_VC = USMC_VC.GetNbSimul();

std::cout << " Berm Case - exec time for control variable MonteCarlo : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_VC).count() << std::endl;
std::cout << " Berm Case - price for control variable MonteCarlo : " << price_US_VC << std::endl;
std::cout << " Berm Case - variance for control variable MonteCarlo : " << var_US_VC << std::endl;
std::cout << " Berm Case - optimal number of simulations for control variable MonteCarlo : " << Optimal_Nb_Simulation_VC << std::endl;

delete ugen_VC;
delete ngen_VC;
delete vectorG_VC;
delete BSamerican_VC;

////////////////// 2.1.2) Control variable (VDC - Laguerre Polynomes n = 2) //////////////////////////////////////

UniformGenerator* vdc_gen_VC = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_VC = new VanDerCorput(11, 1);
Normal* ngen_vdc_VC = new NormalBoxMullerVDC(vdc_gen_VC, vdc_gen_alias_VC, 0., 1.);

GaussianVector* vectorG_vdc_VC = new GaussianVectorCholesky(ngen_vdc_VC, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_VC = new BSEulerND(vectorG_vdc_VC, spot_m, rate);

Bermudean_BasketOption_CV USMC_VC_vdc(N, bsktcall, bsktcallCV, BSamerican_vdc_VC,
	basefunc_Laguerre, wkday, Schedule_exec, price_cf);

clock::time_point start_US_vdc_VC = clock::now();
USMC_VC_vdc.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_VC = clock::now();
clock::duration execution_timeUS_vdc_VC = end_US_vdc_VC - start_US_vdc_VC;
double price_US_vdc_VC = USMC_VC_vdc.GetPrice();
double var_US_vdc_VC = USMC_VC_vdc.GetVariance();

USMC_VC_vdc.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_vdc_VC = USMC_VC_vdc.GetNbSimul();

std::cout << " Berm Case - exec time for control variable MonteCarlo with VDC sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_VC).count() << std::endl;
std::cout << " Berm Case - price for control variable MonteCarlo with VDC sequence : " << price_US_vdc_VC << std::endl;
std::cout << " Berm Case - variance for control variable MonteCarlo  with VDC sequence : " << var_US_vdc_VC << std::endl;
std::cout << " Berm Case - optimal number of simulations for control variable MonteCarlo with VDC sequence : " << Optimal_Nb_Simulation_vdc_VC << std::endl;

delete vdc_gen_VC;
delete vdc_gen_alias_VC;
delete ngen_vdc_VC;
delete vectorG_vdc_VC;
delete BSamerican_vdc_VC;

////////////////// 2.1.3) Control variable (Sobol - Laguerre Polynomes n = 2) //////////////////////////////////////

UniformGenerator* sobol_ugen_CV = new Sobol();
Normal* ngen_sobol_CV = new NormalBoxMuller(sobol_ugen_CV, 0., 1.);

GaussianVector* vectorG_sobol_CV = new GaussianVectorCholesky(ngen_sobol_CV, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_CV = new BSEulerND(vectorG_sobol_CV, spot_m, rate);

Bermudean_BasketOption_CV USMC_sobol_CV(N, bsktcall, bsktcallCV, BSamerican_sobol_CV, basefunc_Laguerre, wkday, Schedule_exec, price_cf);

clock::time_point start_US_sobol_CV = clock::now();
USMC_sobol_CV.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_CV = clock::now();
clock::duration execution_timeUS_sobol_CV = end_US_sobol_CV - start_US_sobol_CV;
double price_US_sobol_CV = USMC_sobol_CV.GetPrice();
double var_US_sobol_CV = USMC_sobol_CV.GetVariance();

USMC_sobol_CV.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_sobol_CV = USMC_sobol_CV.GetNbSimul();

std::cout << " Berm Case - exec time for control variable MonteCarlo with Sobol sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_CV).count() << std::endl;
std::cout << " Berm Case - price for control variable MonteCarlo with Sobol sequence : " << price_US_sobol_CV << std::endl;
std::cout << " Berm Case - variance for control variable MonteCarlo with Sobol sequence : " << var_US_sobol_CV << std::endl;
std::cout << " Berm Case - optimal number of simulations forcontrol variable  MonteCarlo  with Sobol sequence : " << Optimal_Nb_Simulation_sobol_CV << std::endl;

delete sobol_ugen_CV;
delete ngen_sobol_CV;
delete vectorG_sobol_CV;
delete BSamerican_sobol_CV;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// 2.2) Control variable (Ecuyer Combined - Hermite Polynomes n = 2) //////////////////////////////////////


UniformGenerator* ugen_VCH = new EcuyerCombined();
Normal* ngen_VCH = new NormalBoxMuller(ugen_VCH, 0., 1.);

GaussianVector* vectorG_VCH = new GaussianVectorCholesky(ngen_VCH, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_VCH = new BSEulerND(vectorG_VCH, spot_m, rate);

Bermudean_BasketOption_CV  USMC_VCH(N, bsktcall, bsktcallCV, BSamerican_VCH,
	basefunc_Hermite, wkday, Schedule_exec, price_cf);

clock::time_point start_US_VCH = clock::now();
USMC_VCH.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_VCH = clock::now();
clock::duration execution_timeUS_VCH = end_US_VCH - start_US_VCH;
double price_US_VCH = USMC_VCH.GetPrice();
double var_US_VCH = USMC_VCH.GetVariance();

USMC_VCH.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_VCH = USMC_VCH.GetNbSimul();

std::cout << " Berm Case - exec time for control variable MonteCarlo and Hermite polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_VCH).count() << std::endl;
std::cout << " Berm Case - price for control variable MonteCarlo and Hermite polynomes  : " << price_US_VCH << std::endl;
std::cout << " Berm Case - variance for control variable MonteCarlo and Hermite polynomes  : " << var_US_VCH << std::endl;
std::cout << " Berm Case - optimal number of simulations for control variable MonteCarlo and Hermite polynomes  : " << Optimal_Nb_Simulation_VCH << std::endl;

delete ugen_VCH;
delete ngen_VCH;
delete vectorG_VCH;
delete BSamerican_VCH;

////////////////// 2.2.2) Control variable (VDC - Hermite Polynomes n = 2) //////////////////////////////////////

UniformGenerator* vdc_gen_VCH = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_VCH = new VanDerCorput(11, 1);
Normal* ngen_vdc_VCH = new NormalBoxMullerVDC(vdc_gen_VCH, vdc_gen_alias_VCH, 0., 1.);

GaussianVector* vectorG_vdc_VCH = new GaussianVectorCholesky(ngen_vdc_VCH, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_VCH = new BSEulerND(vectorG_vdc_VCH, spot_m, rate);

Bermudean_BasketOption_CV USMC_VC_vdc_Hermite(N, bsktcall, bsktcallCV, BSamerican_vdc_VCH,
	basefunc_Hermite, wkday, Schedule_exec, price_cf);

clock::time_point start_US_vdc_VCH = clock::now();
USMC_VC_vdc_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_VCH = clock::now();
clock::duration execution_timeUS_vdc_VCH = end_US_vdc_VCH - start_US_vdc_VCH;
double price_US_vdc_VCH = USMC_VC_vdc_Hermite.GetPrice();
double var_US_vdc_VCH = USMC_VC_vdc_Hermite.GetVariance();

USMC_VC_vdc_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_vdc_VCH = USMC_VC_vdc_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for control variable MonteCarlo with VDC sequence and Hermite polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_VCH).count() << std::endl;
std::cout << " Berm Case - price for control variable MonteCarlo with VDC sequence and Hermite polynomes  : " << price_US_vdc_VCH << std::endl;
std::cout << " Berm Case - variance for control variable MonteCarlo  with VDC sequence and Hermite polynomes : " << var_US_vdc_VCH << std::endl;
std::cout << " Berm Case - optimal number of simulations for control variable MonteCarlo with VDC sequence and Hermite polynomes : " << Optimal_Nb_Simulation_vdc_VCH << std::endl;

delete vdc_gen_VCH;
delete vdc_gen_alias_VCH;
delete ngen_vdc_VCH;
delete vectorG_vdc_VCH;
delete BSamerican_vdc_VCH;

////////////////// 2.2.3) Control variable (Sobol -Hermite Polynomes n = 2) //////////////////////////////////////

UniformGenerator* sobol_ugen_CVH = new Sobol();
Normal* ngen_sobol_CVH = new NormalBoxMuller(sobol_ugen_CVH, 0., 1.);

GaussianVector* vectorG_sobol_CVH = new GaussianVectorCholesky(ngen_sobol_CVH, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_CVH = new BSEulerND(vectorG_sobol_CVH, spot_m, rate);

Bermudean_BasketOption_CV USMC_sobol_CV_Hermite(N, bsktcall, bsktcallCV, BSamerican_sobol_CVH, basefunc_Hermite, wkday, Schedule_exec, price_cf);

clock::time_point start_US_sobol_CVH = clock::now();
USMC_sobol_CV_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_CVH = clock::now();
clock::duration execution_timeUS_sobol_CVH = end_US_sobol_CVH - start_US_sobol_CVH;
double price_US_sobol_CVH = USMC_sobol_CV_Hermite.GetPrice();
double var_US_sobol_CVH = USMC_sobol_CV_Hermite.GetVariance();

USMC_sobol_CV_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_sobol_CVH = USMC_sobol_CV_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for control variable MonteCarlo with Sobol sequence and Hermite polynomes  : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_CVH).count() << std::endl;
std::cout << " Berm Case - price for control variable MonteCarlo with Sobol sequence and Hermite polynomes  : " << price_US_sobol_CVH << std::endl;
std::cout << " Berm Case - variance for control variable MonteCarlo with Sobol sequence and Hermite polynomes  : " << var_US_sobol_CVH << std::endl;
std::cout << " Berm Case - optimal number of simulations forcontrol variable  MonteCarlo  with Sobol sequence and Hermite polynomes  : " << Optimal_Nb_Simulation_sobol_CVH << std::endl;

delete sobol_ugen_CVH;
delete ngen_sobol_CVH;
delete vectorG_sobol_CVH;
delete BSamerican_sobol_CVH;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// 2.3) Controle variable (Ecuyer Combined - Simples Polynomes n = 2) //////////////////////////////////////


UniformGenerator* ugen_VCS = new EcuyerCombined();
Normal* ngen_VCS = new NormalBoxMuller(ugen_VCS, 0., 1.);

GaussianVector* vectorG_VCS = new GaussianVectorCholesky(ngen_VCS, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_VCS = new BSEulerND(vectorG_VCS, spot_m, rate);

Bermudean_BasketOption_CV USMC_VCS(N, bsktcall, bsktcallCV, BSamerican_VCS,
	basefunc_Simple, wkday, Schedule_exec, price_cf);

clock::time_point start_US_VCS = clock::now();
USMC_VCS.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_VCS = clock::now();
clock::duration execution_timeUS_VCS = end_US_VCS - start_US_VCS;
double price_US_VCS = USMC_VCS.GetPrice();
double var_US_VCS = USMC_VCS.GetVariance();

USMC_VCS.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_VCS = USMC_VCS.GetNbSimul();

std::cout << " Berm Case - exec time for control variable MonteCarlo and Simples polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_VCS).count() << std::endl;
std::cout << " Berm Case - price for control variable MonteCarlo and Simples polynomes  : " << price_US_VCS << std::endl;
std::cout << " Berm Case - variance for control variable MonteCarlo and Simples polynomes  : " << var_US_VCS << std::endl;
std::cout << " Berm Case - optimal number of simulations for control variable MonteCarlo and Simples polynomes  : " << Optimal_Nb_Simulation_VCS << std::endl;

delete ugen_VCS;
delete ngen_VCS;
delete vectorG_VCS;
delete BSamerican_VCS;

////////////////// 2.3.2) Control variable (VDC - Simples Polynomes n = 2) //////////////////////////////////////

UniformGenerator* vdc_gen_VCS = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_VCS = new VanDerCorput(11, 1);
Normal* ngen_vdc_VCS = new NormalBoxMullerVDC(vdc_gen_VCS, vdc_gen_alias_VCS, 0., 1.);

GaussianVector* vectorG_vdc_VCS = new GaussianVectorCholesky(ngen_vdc_VCS, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_VCS = new BSEulerND(vectorG_vdc_VCS, spot_m, rate);

Bermudean_BasketOption_CV USMC_VC_vdc_Simple(N, bsktcall, bsktcallCV, BSamerican_vdc_VCS,
	basefunc_Simple, wkday, Schedule_exec, price_cf);

clock::time_point start_US_vdc_VCS = clock::now();
USMC_VC_vdc_Simple.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_VCS = clock::now();
clock::duration execution_timeUS_vdc_VCS = end_US_vdc_VCS - start_US_vdc_VCS;
double price_US_vdc_VCS = USMC_VC_vdc_Simple.GetPrice();
double var_US_vdc_VCS = USMC_VC_vdc_Simple.GetVariance();

USMC_VC_vdc_Simple.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_vdc_VCS = USMC_VC_vdc_Simple.GetNbSimul();

std::cout << " Berm Case - exec time for control variable MonteCarlo with VDC sequence and Simple polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_VCS).count() << std::endl;
std::cout << " Berm Case - price for control variable MonteCarlo with VDC sequence and Simple polynomes  : " << price_US_vdc_VCS << std::endl;
std::cout << " Berm Case - variance for control variable MonteCarlo  with VDC sequence and Simple polynomes : " << var_US_vdc_VCS << std::endl;
std::cout << " Berm Case - optimal number of simulations for control variable MonteCarlo with VDC sequence and Simple polynomes : " << Optimal_Nb_Simulation_vdc_VCS << std::endl;

delete vdc_gen_VCS;
delete vdc_gen_alias_VCS;
delete ngen_vdc_VCS;
delete vectorG_vdc_VCS;
delete BSamerican_vdc_VCS;

////////////////// 2.3.3) Control variable (Sobol -Simple Polynomes n = 2) //////////////////////////////////////

UniformGenerator* sobol_ugen_CVS = new Sobol();
Normal* ngen_sobol_CVS = new NormalBoxMuller(sobol_ugen_CVS, 0., 1.);

GaussianVector* vectorG_sobol_CVS = new GaussianVectorCholesky(ngen_sobol_CVS, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_CVS = new BSEulerND(vectorG_sobol_CVS, spot_m, rate);

Bermudean_BasketOption_CV USMC_sobol_CV_Simple(N, bsktcall, bsktcallCV, BSamerican_sobol_CVS, basefunc_Simple, wkday, Schedule_exec, price_cf);

clock::time_point start_US_sobol_CVS = clock::now();
USMC_sobol_CV_Simple.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_CVS = clock::now();
clock::duration execution_timeUS_sobol_CVS = end_US_sobol_CVS - start_US_sobol_CVS;
double price_US_sobol_CVS = USMC_sobol_CV_Simple.GetPrice();
double var_US_sobol_CVS = USMC_sobol_CV_Simple.GetVariance();

USMC_sobol_CV_Simple.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_sobol_CVS = USMC_sobol_CV_Simple.GetNbSimul();

std::cout << " Berm Case - exec time for control variable MonteCarlo with Sobol sequence and Simple polynomes  : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_CVS).count() << std::endl;
std::cout << " Berm Case - price for control variable MonteCarlo with Sobol sequence and Simple polynomes  : " << price_US_sobol_CVS << std::endl;
std::cout << " Berm Case - variance for control variable MonteCarlo with Sobol sequence and Simple polynomes  : " << var_US_sobol_CVS << std::endl;
std::cout << " Berm Case - optimal number of simulations forcontrol variable  MonteCarlo  with Sobol sequence and Simple polynomes  : " << Optimal_Nb_Simulation_sobol_CVS << std::endl;

delete sobol_ugen_CVS;
delete ngen_sobol_CVS;
delete vectorG_sobol_CVS;
delete BSamerican_sobol_CVS;

//////////////////////////// Antithetic /////////////////////////////////////////////////////////////////////

//////////////////////// 3.1) Antithetic (Ecuyer Combined - Laguerre polynome n = 2) ////////////////////////////

UniformGenerator* ugen_anti = new EcuyerCombined();
Normal* ngen_anti = new NormalBoxMuller(ugen_anti, 0., 1.);

GaussianVector* vectorG_anti = new GaussianVectorCholesky(ngen_anti, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_anti = new BSEulerNDAntithetic(vectorG_anti, spot_m, rate);

Bermudean_BasketOption_antithetic USMC_anti(N, bsktcall, BSamerican_anti,
	basefunc_Laguerre, wkday, Schedule_exec);

clock::time_point start_US_anti = clock::now();
USMC_anti.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_anti = clock::now();
clock::duration execution_timeUS_anti = end_US_anti - start_US_anti;
double price_US_anti = USMC_anti.GetPrice();
double var_US_anti = USMC_anti.GetVariance();

USMC_anti.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_anti = USMC_anti.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic MonteCarlo : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti).count() << std::endl;
std::cout << " Berm Case - price for antithetic MonteCarlo : " << price_US_anti << std::endl;
std::cout << " Berm Case - variance for antithetic MonteCarlo : " << var_US_anti << std::endl;
std::cout << " Berm Case - optimal number of simulations for antithetic MonteCarlo : " << Optimal_Nb_Simulation_anti << std::endl;

delete ugen_anti;
delete ngen_anti;
delete vectorG_anti;
delete BSamerican_anti;


//////////////////////// 3.1.2) Antithetic (VDC - Laguerre polynome n = 2) ////////////////////////////

UniformGenerator* vdc_gen_anti = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_anti = new VanDerCorput(11, 1);
Normal* ngen_vdc_anti = new NormalBoxMullerVDC(vdc_gen_anti, vdc_gen_alias_anti, 0., 1.);

GaussianVector* vectorG_vdc_anti = new GaussianVectorCholesky(ngen_vdc_anti, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_anti = new BSEulerNDAntithetic(vectorG_vdc_anti, spot_m, rate);

Bermudean_BasketOption_antithetic USMC_vdc_anti(N, bsktcall, BSamerican_vdc_anti,
	basefunc_Laguerre, wkday, Schedule_exec);

clock::time_point start_US_vdc_anti = clock::now();
USMC_vdc_anti.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_anti = clock::now();
clock::duration execution_timeUS_vdc_anti = end_US_vdc_anti - start_US_vdc_anti;
double price_US_vdc_anti = USMC_vdc_anti.GetPrice();
double var_US_vdc_anti = USMC_vdc_anti.GetVariance();

USMC_vdc_anti.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_vdc_anti = USMC_vdc_anti.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic MonteCarlo with VDC sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti).count() << std::endl;
std::cout << " Berm Case - price for antithetic MonteCarlo with VDC sequence : " << price_US_vdc_anti << std::endl;
std::cout << " Berm Case - variance for antithetic MonteCarlo with VDC sequence : " << var_US_vdc_anti << std::endl;
std::cout << " Berm Case - optimal number of simulations for antithetic MonteCarlo with VDC sequence : " << Optimal_Nb_Simulation_vdc_anti << std::endl;

delete vdc_gen_anti;
delete vdc_gen_alias_anti;
delete ngen_vdc_anti;
delete vectorG_vdc_anti;
delete BSamerican_vdc_anti;

//////////////////////// 3.1.3) Antithetic (Sobol - Laguerre polynome n = 2) ////////////////////////////

UniformGenerator* gen_sobol_anti = new Sobol();
Normal* ngen_sobol_anti = new NormalBoxMuller(gen_sobol_anti, 0., 1.);

GaussianVector* vectorG_sobol_anti = new GaussianVectorCholesky(ngen_sobol_anti, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_anti = new BSEulerNDAntithetic(vectorG_sobol_anti, spot_m, rate);

Bermudean_BasketOption_antithetic USMC_sobol_anti(N, bsktcall, BSamerican_sobol_anti,
	basefunc_Laguerre, wkday, Schedule_exec);

clock::time_point start_US_sobol_anti = clock::now();
USMC_sobol_anti.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_anti = clock::now();
clock::duration execution_timeUS_sobol_anti = end_US_sobol_anti - start_US_sobol_anti;
double price_US_sobol_anti = USMC_sobol_anti.GetPrice();
double var_US_sobol_anti = USMC_sobol_anti.GetVariance();

USMC_sobol_anti.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_sobol_anti = USMC_sobol_anti.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic MonteCarlo with Sobol sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti).count() << std::endl;
std::cout << " Berm Case - price for antithetic MonteCarlo with Sobol sequence : " << price_US_sobol_anti << std::endl;
std::cout << " Berm Case - variance for antithetic MonteCarlo with Sobol sequence : " << var_US_sobol_anti << std::endl;
std::cout << " Berm Case - optimal number of simulations for antithetic MonteCarlo with Sobol sequence : " << Optimal_Nb_Simulation_sobol_anti << std::endl;

delete gen_sobol_anti;
delete ngen_sobol_anti;
delete vectorG_sobol_anti;
delete BSamerican_sobol_anti;

//////////////////////// 3.2) Antithetic (Ecuyer Combined - Hermite polynome n = 2) ////////////////////////////

UniformGenerator* ugen_antiH = new EcuyerCombined();
Normal* ngen_antiH = new NormalBoxMuller(ugen_antiH, 0., 1.);

GaussianVector* vectorG_antiH = new GaussianVectorCholesky(ngen_antiH, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_antiH = new BSEulerNDAntithetic(vectorG_antiH, spot_m, rate);

Bermudean_BasketOption_antithetic USMC_anti_Hermite(N, bsktcall, BSamerican_antiH,
	basefunc_Hermite, wkday, Schedule_exec);

clock::time_point start_US_antiH = clock::now();
USMC_anti_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_antiH = clock::now();
clock::duration execution_timeUS_antiH = end_US_antiH - start_US_antiH;
double price_US_antiH = USMC_anti_Hermite.GetPrice();
double var_US_antiH = USMC_anti_Hermite.GetVariance();

USMC_anti_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_antiH = USMC_anti_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic MonteCarlo with Hermite polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_antiH).count() << std::endl;
std::cout << " Berm Case - price for antithetic MonteCarlo with Hermite polynome : " << price_US_antiH << std::endl;
std::cout << " Berm Case - variance for antithetic MonteCarlo with Hermite polynome : " << var_US_antiH << std::endl;
std::cout << " Berm Case - optimal number of simulations for antithetic MonteCarlo with Hermite polynome : " << Optimal_Nb_Simulation_antiH << std::endl;

delete ugen_antiH;
delete ngen_antiH;
delete vectorG_antiH;
delete BSamerican_antiH;


//////////////////////// 3.2.2) Antithetic (VDC - Hermite polynome n = 2) ////////////////////////////

UniformGenerator* vdc_gen_antiH = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_antiH = new VanDerCorput(11, 1);
Normal* ngen_vdc_antiH = new NormalBoxMullerVDC(vdc_gen_antiH, vdc_gen_alias_antiH, 0., 1.);

GaussianVector* vectorG_vdc_antiH = new GaussianVectorCholesky(ngen_vdc_antiH, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_antiH = new BSEulerNDAntithetic(vectorG_vdc_antiH, spot_m, rate);

Bermudean_BasketOption_antithetic USMC_vdc_anti_Hermite(N, bsktcall, BSamerican_vdc_antiH,
	basefunc_Hermite, wkday, Schedule_exec);

clock::time_point start_US_vdc_antiH = clock::now();
USMC_vdc_anti_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_antiH = clock::now();
clock::duration execution_timeUS_vdc_antiH = end_US_vdc_antiH - start_US_vdc_antiH;
double price_US_vdc_antiH = USMC_vdc_anti_Hermite.GetPrice();
double var_US_vdc_antiH = USMC_vdc_anti_Hermite.GetVariance();

USMC_vdc_anti_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_vdc_antiH = USMC_vdc_anti_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic MonteCarlo with VDC sequence and Hermite Polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_antiH).count() << std::endl;
std::cout << " Berm Case - price for antithetic MonteCarlo with VDC sequence and Hermite Polynome : " << price_US_vdc_antiH << std::endl;
std::cout << " Berm Case - variance for antithetic MonteCarlo with VDC sequence and Hermite Polynome : " << var_US_vdc_antiH << std::endl;
std::cout << " Berm Case - optimal number of simulations for antithetic MonteCarlo with VDC sequence and Hermite Polynome : " << Optimal_Nb_Simulation_vdc_antiH << std::endl;

delete vdc_gen_antiH;
delete vdc_gen_alias_antiH;
delete ngen_vdc_antiH;
delete vectorG_vdc_antiH;
delete BSamerican_vdc_antiH;

//////////////////////// 3.2.3) Antithetic (Sobol - Hermite polynome n = 2) ////////////////////////////

UniformGenerator* gen_sobol_antiH = new Sobol();
Normal* ngen_sobol_antiH = new NormalBoxMuller(gen_sobol_antiH, 0., 1.);

GaussianVector* vectorG_sobol_antiH = new GaussianVectorCholesky(ngen_sobol_antiH, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_antiH = new BSEulerNDAntithetic(vectorG_sobol_antiH, spot_m, rate);

Bermudean_BasketOption_antithetic USMC_sobol_anti_Hermite(N, bsktcall, BSamerican_sobol_antiH,
	basefunc_Hermite, wkday, Schedule_exec);

clock::time_point start_US_sobol_antiH = clock::now();
USMC_sobol_anti_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_antiH = clock::now();
clock::duration execution_timeUS_sobol_antiH = end_US_sobol_antiH - start_US_sobol_antiH;
double price_US_sobol_antiH = USMC_sobol_anti_Hermite.GetPrice();
double var_US_sobol_antiH = USMC_sobol_anti_Hermite.GetVariance();

USMC_sobol_anti_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_sobol_antiH = USMC_sobol_anti_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic MonteCarlo with Sobol sequence and Hermite Polynome: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_antiH).count() << std::endl;
std::cout << " Berm Case - price for antithetic MonteCarlo with Sobol sequence and Hermite Polynome : " << price_US_sobol_antiH << std::endl;
std::cout << " Berm Case - variance for antithetic MonteCarlo with Sobol sequence and Hermite Polynome : " << var_US_sobol_antiH << std::endl;
std::cout << " Berm Case - optimal number of simulations for antithetic MonteCarlo with Sobol sequence and Hermite Polynome : " << Optimal_Nb_Simulation_sobol_antiH << std::endl;

delete gen_sobol_antiH;
delete ngen_sobol_antiH;
delete vectorG_sobol_antiH;
delete BSamerican_sobol_antiH;

//////////////////////// 3.3) Antithetic (Ecuyer Combined - Simple polynome n = 2) ////////////////////////////

UniformGenerator* ugen_antiS = new EcuyerCombined();
Normal* ngen_antiS = new NormalBoxMuller(ugen_antiS, 0., 1.);

GaussianVector* vectorG_antiS = new GaussianVectorCholesky(ngen_antiS, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_antiS = new BSEulerNDAntithetic(vectorG_antiS, spot_m, rate);

Bermudean_BasketOption_antithetic USMC_anti_Simple(N, bsktcall, BSamerican_antiS,
	basefunc_Simple, wkday, Schedule_exec);

clock::time_point start_US_antiS = clock::now();
USMC_anti_Simple.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_antiS = clock::now();
clock::duration execution_timeUS_antiS = end_US_antiS - start_US_antiS;
double price_US_antiS = USMC_anti_Simple.GetPrice();
double var_US_antiS = USMC_anti_Simple.GetVariance();

USMC_anti_Simple.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_antiS = USMC_anti_Simple.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic MonteCarlo with Simple polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_antiS).count() << std::endl;
std::cout << " Berm Case - price for antithetic MonteCarlo with Simple polynome : " << price_US_antiS << std::endl;
std::cout << " Berm Case - variance for antithetic MonteCarlo with Simple polynome : " << var_US_antiS << std::endl;
std::cout << " Berm Case - optimal number of simulations for antithetic MonteCarlo with Simple polynome : " << Optimal_Nb_Simulation_antiS << std::endl;

delete ugen_antiS;
delete ngen_antiS;
delete vectorG_antiS;
delete BSamerican_antiS;


//////////////////////// 3.3.2) Antithetic (VDC - Simple polynome n = 2) ////////////////////////////

UniformGenerator* vdc_gen_antiS = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_antiS = new VanDerCorput(11, 1);
Normal* ngen_vdc_antiS = new NormalBoxMullerVDC(vdc_gen_antiS, vdc_gen_alias_antiS, 0., 1.);

GaussianVector* vectorG_vdc_antiS = new GaussianVectorCholesky(ngen_vdc_antiS, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_antiS = new BSEulerNDAntithetic(vectorG_vdc_antiS, spot_m, rate);

Bermudean_BasketOption_antithetic USMC_vdc_anti_Simple(N, bsktcall, BSamerican_vdc_antiS,
	basefunc_Simple, wkday, Schedule_exec);

clock::time_point start_US_vdc_antiS = clock::now();
USMC_vdc_anti_Simple.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_antiS = clock::now();
clock::duration execution_timeUS_vdc_antiS = end_US_vdc_antiS - start_US_vdc_antiS;
double price_US_vdc_antiS = USMC_vdc_anti_Simple.GetPrice();
double var_US_vdc_antiS = USMC_vdc_anti_Simple.GetVariance();

USMC_vdc_anti_Simple.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_vdc_antiS = USMC_vdc_anti_Simple.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic MonteCarlo with VDC sequence and Simple Polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_antiS).count() << std::endl;
std::cout << " Berm Case - price for antithetic MonteCarlo with VDC sequence and Simple Polynome : " << price_US_vdc_antiS << std::endl;
std::cout << " Berm Case - variance for antithetic MonteCarlo with VDC sequence and Simple Polynome : " << var_US_vdc_antiS << std::endl;
std::cout << " Berm Case - optimal number of simulations for antithetic MonteCarlo with VDC sequence and Simple Polynome : " << Optimal_Nb_Simulation_vdc_antiS << std::endl;

delete vdc_gen_antiS;
delete vdc_gen_alias_antiS;
delete ngen_vdc_antiS;
delete vectorG_vdc_antiS;
delete BSamerican_vdc_antiS;

//////////////////////// 3.3.3) Antithetic (Sobol - Simple polynome n = 2) ////////////////////////////

UniformGenerator* gen_sobol_antiS = new Sobol();
Normal* ngen_sobol_antiS = new NormalBoxMuller(gen_sobol_antiS, 0., 1.);

GaussianVector* vectorG_sobol_antiS = new GaussianVectorCholesky(ngen_sobol_antiS, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_antiS = new BSEulerNDAntithetic(vectorG_sobol_antiS, spot_m, rate);

Bermudean_BasketOption_antithetic USMC_sobol_anti_Simple(N, bsktcall, BSamerican_sobol_antiS,
	basefunc_Simple, wkday, Schedule_exec);

clock::time_point start_US_sobol_antiS = clock::now();
USMC_sobol_anti_Simple.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_antiS = clock::now();
clock::duration execution_timeUS_sobol_antiS = end_US_sobol_antiS - start_US_sobol_antiS;
double price_US_sobol_antiS = USMC_sobol_anti_Simple.GetPrice();
double var_US_sobol_antiS = USMC_sobol_anti_Simple.GetVariance();

USMC_sobol_anti_Simple.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_sobol_antiS = USMC_sobol_anti_Simple.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic MonteCarlo with Sobol sequence and Simple Polynome: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_antiS).count() << std::endl;
std::cout << " Berm Case - price for antithetic MonteCarlo with Sobol sequence and Simple Polynome : " << price_US_sobol_antiS << std::endl;
std::cout << " Berm Case - variance for antithetic MonteCarlo with Sobol sequence and Simple Polynome : " << var_US_sobol_antiS << std::endl;
std::cout << " Berm Case - optimal number of simulations for antithetic MonteCarlo with Sobol sequence and Simple Polynome : " << Optimal_Nb_Simulation_sobol_antiS << std::endl;

delete gen_sobol_antiS;
delete ngen_sobol_antiS;
delete vectorG_sobol_antiS;
delete BSamerican_sobol_antiS;

//////////////////////////// VC + Antithetic /////////////////////////////////////////////////////////////////////

//////////////////////// 4.1) VC + Antithetic (Ecuyer Combined - Laguerre polynome n = 2) ////////////////////////////

UniformGenerator* ugen_anti_cv = new EcuyerCombined();
Normal* ngen_anti_cv = new NormalBoxMuller(ugen_anti_cv, 0., 1.);

GaussianVector* vectorG_anti_cv = new GaussianVectorCholesky(ngen_anti_cv, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_anti_cv = new BSEulerNDAntithetic(vectorG_anti_cv, spot_m, rate);

Bermudean_BasketOption_antithetic_CV USMC_anti_CV(N, bsktcall, bsktcallCV, BSamerican_anti_cv, basefunc_Laguerre, wkday, Schedule_exec, price_cf);


clock::time_point start_US_anti_cv = clock::now();
USMC_anti_CV.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_anti_cv = clock::now();
clock::duration execution_timeUS_anti_cv = end_US_anti_cv - start_US_anti_cv;
double price_US_anti_cv = USMC_anti_CV.GetPrice();
double var_US_anti_cv = USMC_anti_CV.GetVariance();

USMC_anti_CV.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_anti_cv = USMC_anti_CV.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic and control variable MonteCarlo : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti_cv).count() << std::endl;
std::cout << " Berm Case - price for antithetic and control variable MonteCarlo : " << price_US_anti_cv << std::endl;
std::cout << " Berm Case - variance for antithetic and control variable MonteCarlo : " << var_US_anti_cv << std::endl;
std::cout << " Berm Case - optimal number of simulations forantithetic and control variable MonteCarlo : " << Optimal_Nb_Simulation_anti_cv << std::endl;

delete ugen_anti_cv;
delete ngen_anti_cv;
delete vectorG_anti_cv;
delete BSamerican_anti_cv;

//////////////////////// 4.1.2) VC + Antithetic (VDC - Laguerre polynome n = 2) ////////////////////////////

UniformGenerator* vdc_gen_anti_cv = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_anti_cv = new VanDerCorput(11, 1);
Normal* ngen_vdc_anti_cv = new NormalBoxMullerVDC(vdc_gen_anti_cv, vdc_gen_alias_anti_cv, 0., 1.);

GaussianVector* vectorG_vdc_anti_cv = new GaussianVectorCholesky(ngen_vdc_anti_cv, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_anti_cv = new BSEulerNDAntithetic(vectorG_vdc_anti_cv, spot_m, rate);

Bermudean_BasketOption_antithetic_CV USMC_vdc_anti_CV(N, bsktcall, bsktcallCV, BSamerican_vdc_anti_cv, basefunc_Laguerre, wkday, Schedule_exec, price_cf);


clock::time_point start_US_vdc_anti_cv = clock::now();
USMC_vdc_anti_CV.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_anti_cv = clock::now();
clock::duration execution_timeUS_vdc_anti_cv = end_US_vdc_anti_cv - start_US_vdc_anti_cv;
double price_US_vdc_anti_cv = USMC_vdc_anti_CV.GetPrice();
double var_US_vdc_anti_cv = USMC_vdc_anti_CV.GetVariance();

USMC_vdc_anti_CV.OptimalNbSimul(tolerated_error);
size_t Optimal_vdc_Nb_Simulation_anti_cv = USMC_vdc_anti_CV.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic and control variable MonteCarlo with VDC sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti_cv).count() << std::endl;
std::cout << " Berm Case - price for antithetic and control variable MonteCarlo  with VDC sequence: " << price_US_vdc_anti_cv << std::endl;
std::cout << " Berm Case - variance for antithetic and control variable MonteCarlo with VDC sequence : " << var_US_vdc_anti_cv << std::endl;
std::cout << " Berm Case - optimal number of simulations forantithetic and control variable MonteCarlo with VDC sequence : " << Optimal_vdc_Nb_Simulation_anti_cv << std::endl;

delete vdc_gen_anti_cv;
delete vdc_gen_alias_anti_cv;
delete ngen_vdc_anti_cv;
delete vectorG_vdc_anti_cv;
delete BSamerican_vdc_anti_cv;

//////////////////////// 4.1.2) VC + Antithetic (Sobol - Laguerre polynome n = 2) ////////////////////////////

UniformGenerator* sobol_gen_anti_cv = new Sobol();
Normal* ngen_sobol_anti_cv = new NormalBoxMuller(sobol_gen_anti_cv, 0., 1.);

GaussianVector* vectorG_sobol_anti_cv = new GaussianVectorCholesky(ngen_sobol_anti_cv, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_anti_cv = new BSEulerNDAntithetic(vectorG_sobol_anti_cv, spot_m, rate);

Bermudean_BasketOption_antithetic_CV USMC_sobol_anti_CV(N, bsktcall, bsktcallCV, BSamerican_sobol_anti_cv, basefunc_Laguerre, wkday, Schedule_exec, price_cf);


clock::time_point start_US_sobol_anti_cv = clock::now();
USMC_sobol_anti_CV.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_anti_cv = clock::now();
clock::duration execution_timeUS_sobol_anti_cv = end_US_sobol_anti_cv - start_US_sobol_anti_cv;
double price_US_sobol_anti_cv = USMC_sobol_anti_CV.GetPrice();
double var_US_sobol_anti_cv = USMC_sobol_anti_CV.GetVariance();

USMC_sobol_anti_CV.OptimalNbSimul(tolerated_error);
size_t Optimal_sobol_Nb_Simulation_anti_cv = USMC_sobol_anti_CV.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic and control variable MonteCarlo with Sobol sequence : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti_cv).count() << std::endl;
std::cout << " Berm Case - price for antithetic and control variable MonteCarlo  with Sobol sequence: " << price_US_sobol_anti_cv << std::endl;
std::cout << " Berm Case - variance for antithetic and control variable MonteCarlo with Sobol sequence : " << var_US_sobol_anti_cv << std::endl;
std::cout << " Berm Case - optimal number of simulations forantithetic and control variable MonteCarlo with Sobol sequence : " << Optimal_sobol_Nb_Simulation_anti_cv << std::endl;

delete sobol_gen_anti_cv;
delete ngen_sobol_anti_cv;
delete vectorG_sobol_anti_cv;
delete BSamerican_sobol_anti_cv;

//////////////////////// 4.2) VC + Antithetic (Ecuyer Combined - Hermite polynome n = 2) ////////////////////////////

UniformGenerator* ugen_anti_cvH = new EcuyerCombined();
Normal* ngen_anti_cvH = new NormalBoxMuller(ugen_anti_cvH, 0., 1.);

GaussianVector* vectorG_anti_cvH = new GaussianVectorCholesky(ngen_anti_cvH, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_anti_cvH = new BSEulerNDAntithetic(vectorG_anti_cvH, spot_m, rate);

Bermudean_BasketOption_antithetic_CV USMC_anti_CV_Hermite(N, bsktcall, bsktcallCV, BSamerican_anti_cvH, basefunc_Hermite, wkday, Schedule_exec, price_cf);


clock::time_point start_US_anti_cvH = clock::now();
USMC_anti_CV_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_anti_cvH = clock::now();
clock::duration execution_timeUS_anti_cvH = end_US_anti_cvH - start_US_anti_cvH;
double price_US_anti_cvH = USMC_anti_CV_Hermite.GetPrice();
double var_US_anti_cvH = USMC_anti_CV_Hermite.GetVariance();

USMC_anti_CV_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_anti_cvH = USMC_anti_CV_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic and control variable MonteCarlo and Hermite polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti_cvH).count() << std::endl;
std::cout << " Berm Case - price for antithetic and control variable MonteCarlo and Hermite polynome : " << price_US_anti_cvH << std::endl;
std::cout << " Berm Case - variance for antithetic and control variable MonteCarlo and Hermite polynome : " << var_US_anti_cvH << std::endl;
std::cout << " Berm Case - optimal number of simulations forantithetic and control variable MonteCarlo and Hermite Polynome : " << Optimal_Nb_Simulation_anti_cvH << std::endl;

delete ugen_anti_cvH;
delete ngen_anti_cvH;
delete vectorG_anti_cvH;
delete BSamerican_anti_cvH;

//////////////////////// 4.2.2) VC + Antithetic (VDC - Hermite polynome n = 2) ////////////////////////////

UniformGenerator* vdc_gen_anti_cvH = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_anti_cvH = new VanDerCorput(11, 1);
Normal* ngen_vdc_anti_cvH = new NormalBoxMullerVDC(vdc_gen_anti_cvH, vdc_gen_alias_anti_cvH, 0., 1.);

GaussianVector* vectorG_vdc_anti_cvH = new GaussianVectorCholesky(ngen_vdc_anti_cvH, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_anti_cvH = new BSEulerNDAntithetic(vectorG_vdc_anti_cvH, spot_m, rate);

Bermudean_BasketOption_antithetic_CV USMC_vdc_anti_CV_Hermite(N, bsktcall, bsktcallCV, BSamerican_vdc_anti_cvH, basefunc_Hermite, wkday, Schedule_exec, price_cf);


clock::time_point start_US_vdc_anti_cvH = clock::now();
USMC_vdc_anti_CV_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_anti_cvH = clock::now();
clock::duration execution_timeUS_vdc_anti_cvH = end_US_vdc_anti_cvH - start_US_vdc_anti_cvH;
double price_US_vdc_anti_cvH = USMC_vdc_anti_CV_Hermite.GetPrice();
double var_US_vdc_anti_cvH = USMC_vdc_anti_CV_Hermite.GetVariance();

USMC_vdc_anti_CV_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_vdc_Nb_Simulation_anti_cvH = USMC_vdc_anti_CV_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic and control variable MonteCarlo with VDC sequence and Hermite polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti_cvH).count() << std::endl;
std::cout << " Berm Case - price for antithetic and control variable MonteCarlo  with VDC sequence and Hermite polynome : " << price_US_vdc_anti_cvH << std::endl;
std::cout << " Berm Case - variance for antithetic and control variable MonteCarlo with VDC sequence and Hermite polynome : " << var_US_vdc_anti_cvH << std::endl;
std::cout << " Berm Case - optimal number of simulations forantithetic and control variable MonteCarlo with VDC sequence and Hermite polynome : " << Optimal_vdc_Nb_Simulation_anti_cvH << std::endl;

delete vdc_gen_anti_cvH;
delete vdc_gen_alias_anti_cvH;
delete ngen_vdc_anti_cvH;
delete vectorG_vdc_anti_cvH;
delete BSamerican_vdc_anti_cvH;

//////////////////////// 4.2.3) VC + Antithetic (Sobol - Hermite polynome n = 2) ////////////////////////////

UniformGenerator* sobol_gen_anti_cvH = new Sobol();
Normal* ngen_sobol_anti_cvH = new NormalBoxMuller(sobol_gen_anti_cvH, 0., 1.);

GaussianVector* vectorG_sobol_anti_cvH = new GaussianVectorCholesky(ngen_sobol_anti_cvH, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_anti_cvH = new BSEulerNDAntithetic(vectorG_sobol_anti_cvH, spot_m, rate);

Bermudean_BasketOption_antithetic_CV USMC_sobol_anti_CV_Hermite(N, bsktcall, bsktcallCV, BSamerican_sobol_anti_cvH, basefunc_Hermite, wkday, Schedule_exec, price_cf);


clock::time_point start_US_sobol_anti_cvH = clock::now();
USMC_sobol_anti_CV_Hermite.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_anti_cvH = clock::now();
clock::duration execution_timeUS_sobol_anti_cvH = end_US_sobol_anti_cvH - start_US_sobol_anti_cvH;
double price_US_sobol_anti_cvH = USMC_sobol_anti_CV_Hermite.GetPrice();
double var_US_sobol_anti_cvH = USMC_sobol_anti_CV_Hermite.GetVariance();

USMC_sobol_anti_CV_Hermite.OptimalNbSimul(tolerated_error);
size_t Optimal_sobol_Nb_Simulation_anti_cvH = USMC_sobol_anti_CV_Hermite.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic and control variable MonteCarlo with Sobol sequence and Hermite polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti_cvH).count() << std::endl;
std::cout << " Berm Case - price for antithetic and control variable MonteCarlo  with Sobol sequence and Hermite polynome : " << price_US_sobol_anti_cvH << std::endl;
std::cout << " Berm Case - variance for antithetic and control variable MonteCarlo with Sobol sequence and Hermite polynome : " << var_US_sobol_anti_cvH << std::endl;
std::cout << " Berm Case - optimal number of simulations forantithetic and control variable MonteCarlo with Sobol sequence and Hermite polynome : " << Optimal_sobol_Nb_Simulation_anti_cvH << std::endl;

delete sobol_gen_anti_cvH;
delete ngen_sobol_anti_cvH;
delete vectorG_sobol_anti_cvH;
delete BSamerican_sobol_anti_cvH;

//////////////////////// 4.3) VC + Antithetic (Ecuyer Combined - simple polynome n = 2) ////////////////////////////

UniformGenerator* ugen_anti_cvS = new EcuyerCombined();
Normal* ngen_anti_cvS = new NormalBoxMuller(ugen_anti_cvS, 0., 1.);

GaussianVector* vectorG_anti_cvS = new GaussianVectorCholesky(ngen_anti_cvS, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_anti_cvS = new BSEulerNDAntithetic(vectorG_anti_cvS, spot_m, rate);

Bermudean_BasketOption_antithetic_CV USMC_anti_CV_Simple(N, bsktcall, bsktcallCV, BSamerican_anti_cvS, basefunc_Simple, wkday, Schedule_exec, price_cf);


clock::time_point start_US_anti_cvS = clock::now();
USMC_anti_CV_Simple.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_anti_cvS = clock::now();
clock::duration execution_timeUS_anti_cvS = end_US_anti_cvS - start_US_anti_cvS;
double price_US_anti_cvS = USMC_anti_CV_Simple.GetPrice();
double var_US_anti_cvS = USMC_anti_CV_Simple.GetVariance();

USMC_anti_CV_Simple.OptimalNbSimul(tolerated_error);
size_t Optimal_Nb_Simulation_anti_cvS = USMC_anti_CV_Simple.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic and control variable MonteCarlo and Simple polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti_cvS).count() << std::endl;
std::cout << " Berm Case - price for antithetic and control variable MonteCarlo and Simple polynome : " << price_US_anti_cvS << std::endl;
std::cout << " Berm Case - variance for antithetic and control variable MonteCarlo and Simple polynome : " << var_US_anti_cvS << std::endl;
std::cout << " Berm Case - optimal number of simulations forantithetic and control variable MonteCarlo and Simple Polynome : " << Optimal_Nb_Simulation_anti_cvS << std::endl;

delete ugen_anti_cvS;
delete ngen_anti_cvS;
delete vectorG_anti_cvS;
delete BSamerican_anti_cvS;

//////////////////////// 4.3.2) VC + Antithetic (VDC - Simple polynome n = 2) ////////////////////////////

UniformGenerator* vdc_gen_anti_cvS = new VanDerCorput(2, 1);
UniformGenerator* vdc_gen_alias_anti_cvS = new VanDerCorput(11, 1);
Normal* ngen_vdc_anti_cvS = new NormalBoxMullerVDC(vdc_gen_anti_cvS, vdc_gen_alias_anti_cvS, 0., 1.);

GaussianVector* vectorG_vdc_anti_cvS = new GaussianVectorCholesky(ngen_vdc_anti_cvS, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_vdc_anti_cvS = new BSEulerNDAntithetic(vectorG_vdc_anti_cvS, spot_m, rate);

Bermudean_BasketOption_antithetic_CV USMC_vdc_anti_CV_Simple(N, bsktcall, bsktcallCV, BSamerican_vdc_anti_cvS, basefunc_Simple, wkday, Schedule_exec, price_cf);


clock::time_point start_US_vdc_anti_cvS = clock::now();
USMC_vdc_anti_CV_Simple.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_vdc_anti_cvS = clock::now();
clock::duration execution_timeUS_vdc_anti_cvS = end_US_vdc_anti_cvS - start_US_vdc_anti_cvS;
double price_US_vdc_anti_cvS = USMC_vdc_anti_CV_Simple.GetPrice();
double var_US_vdc_anti_cvS = USMC_vdc_anti_CV_Simple.GetVariance();

USMC_vdc_anti_CV_Simple.OptimalNbSimul(tolerated_error);
size_t Optimal_vdc_Nb_Simulation_anti_cvS = USMC_vdc_anti_CV_Simple.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic and control variable MonteCarlo with VDC sequence and simple polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti_cvS).count() << std::endl;
std::cout << " Berm Case - price for antithetic and control variable MonteCarlo  with VDC sequence and simple polynome : " << price_US_vdc_anti_cvS << std::endl;
std::cout << " Berm Case - variance for antithetic and control variable MonteCarlo with VDC sequence and simple polynome : " << var_US_vdc_anti_cvS << std::endl;
std::cout << " Berm Case - optimal number of simulations forantithetic and control variable MonteCarlo with VDC sequence and simple polynome : " << Optimal_vdc_Nb_Simulation_anti_cvS << std::endl;

delete vdc_gen_anti_cvS;
delete vdc_gen_alias_anti_cvS;
delete ngen_vdc_anti_cvS;
delete vectorG_vdc_anti_cvS;
delete BSamerican_vdc_anti_cvS;

//////////////////////// 4.3.3) VC + Antithetic (Sobol - simple polynome n = 2) ////////////////////////////

UniformGenerator* sobol_gen_anti_cvS = new Sobol();
Normal* ngen_sobol_anti_cvS = new NormalBoxMuller(sobol_gen_anti_cvS, 0., 1.);

GaussianVector* vectorG_sobol_anti_cvS = new GaussianVectorCholesky(ngen_sobol_anti_cvS, Sigma, Correl, CovarMatrix);
RandomProcess* BSamerican_sobol_anti_cvS = new BSEulerNDAntithetic(vectorG_sobol_anti_cvS, spot_m, rate);

Bermudean_BasketOption_antithetic_CV USMC_sobol_anti_CV_Simple(N, bsktcall, bsktcallCV, BSamerican_sobol_anti_cvS, basefunc_Simple, wkday, Schedule_exec, price_cf);


clock::time_point start_US_sobol_anti_cvS = clock::now();
USMC_sobol_anti_CV_Simple.Simulate(startTime, endTime, nbsteps);
clock::time_point end_US_sobol_anti_cvS = clock::now();
clock::duration execution_timeUS_sobol_anti_cvS = end_US_sobol_anti_cvS - start_US_sobol_anti_cvS;
double price_US_sobol_anti_cvS = USMC_sobol_anti_CV_Simple.GetPrice();
double var_US_sobol_anti_cvS = USMC_sobol_anti_CV_Simple.GetVariance();

USMC_sobol_anti_CV_Simple.OptimalNbSimul(tolerated_error);
size_t Optimal_sobol_Nb_Simulation_anti_cvS = USMC_sobol_anti_CV_Simple.GetNbSimul();

std::cout << " Berm Case - exec time for antithetic and control variable MonteCarlo with Sobol sequence and simple polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti_cvS).count() << std::endl;
std::cout << " Berm Case - price for antithetic and control variable MonteCarlo  with Sobol sequence and simple polynome : " << price_US_sobol_anti_cvS << std::endl;
std::cout << " Berm Case - variance for antithetic and control variable MonteCarlo with Sobol sequence and Simple polynome : " << var_US_sobol_anti_cvS << std::endl;
std::cout << " Berm Case - optimal number of simulations forantithetic and control variable MonteCarlo with Sobol sequence and Simple polynome : " << Optimal_sobol_Nb_Simulation_anti_cvS << std::endl;

delete sobol_gen_anti_cvS;
delete ngen_sobol_anti_cvS;
delete vectorG_sobol_anti_cvS;
delete BSamerican_sobol_anti_cvS; */



/////////////////////////::::test for higher orders of polynomes::::::////////////////////////////////////////////////////////



//// need to choose the polynome basis with which one will approxime the continuation value of the option
//	std::vector<basis_functions*> basefunc_Laguerre_6;
//
//	basefunc_Laguerre_6.push_back(new Poly_Laguerre(0));
//	basefunc_Laguerre_6.push_back(new Poly_Laguerre(1));
//	basefunc_Laguerre_6.push_back(new Poly_Laguerre(2));
//	basefunc_Laguerre_6.push_back(new Poly_Laguerre(3));
//	basefunc_Laguerre_6.push_back(new Poly_Laguerre(4));
//	basefunc_Laguerre_6.push_back(new Poly_Laguerre(5));
//	basefunc_Laguerre_6.push_back(new Poly_Laguerre(6));
//
//	//Create the payoff of the options 
//
//	PayOffBasket* bsktcall = new PayOffBasketCall(W, spot_m, 100.);
//
//	/////////////////// NO VARIANCE REDUCTION /////////////////////////////
//
//
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	////////////////// 1.4) Vanille (Ecuyer Combined - Laguerre Polynomes n = 6) //////////////////////////////////////
//
//	UniformGenerator* ugen10 = new EcuyerCombined();
//	Normal* ngen10 = new NormalBoxMuller(ugen10, 0., 1.);
//
//	GaussianVector* vectorG10 = new GaussianVectorCholesky(ngen10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican10 = new BSEulerND(vectorG10, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_10(N, bsktcall, BSamerican10, basefunc_Laguerre_6);
//
//	clock::time_point start_US_10 = clock::now();
//	USMC_10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_10 = clock::now();
//	clock::duration execution_timeUS_10 = end_US_10 - start_US_10;
//	double price_US_10 = USMC_10.GetPrice();
//	double var_US_10 = USMC_10.GetVariance();
//
//	USMC_10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_10 = USMC_10.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with 6 order of Laguerre polynomes: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_10).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with 6 order of Laguerre polynomes : " << price_US_10 << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo with 6 order of Laguerre polynomes : " << var_US_10 << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo with 6 order of Laguerre polynomes : " << Optimal_Nb_Simulation_10 << std::endl;
//
//	delete ugen10;
//	delete ngen10;
//	delete vectorG10;
//	delete BSamerican10;
//
//	////////////////// 1.4.2) Vanille (VDC - Laguerre Polynomes n = 6) //////////////////////////////////////
//
//	UniformGenerator* vdc_gen_10 = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_10 = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_10 = new NormalBoxMullerVDC(vdc_gen_10, vdc_gen_alias_10, 0., 1.);
//
//	GaussianVector* vectorG_vdc_10 = new GaussianVectorCholesky(ngen_vdc_10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_10 = new BSEulerND(vectorG_vdc_10, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_vdc_10(N, bsktcall, BSamerican_vdc_10, basefunc_Laguerre_6);
//
//	clock::time_point start_US_vdc_10 = clock::now();
//	USMC_vdc_10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_10 = clock::now();
//	clock::duration execution_timeUS_vdc_10 = end_US_vdc_10 - start_US_vdc_10;
//	double price_US_vdc_10 = USMC_vdc_10.GetPrice();
//	double var_US_vdc_10 = USMC_vdc_10.GetVariance();
//
//	USMC_vdc_10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_10 = USMC_vdc_10.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with VDC sequence and 6 order of Laguerre polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_10).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with VDC sequence and 6 order of Laguerre polynomes : " << price_US_vdc_10 << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo  with VDC sequence and 6 order of Laguerre polynomes : " << var_US_vdc_10 << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo  with VDC sequence and 6 order of Laguerre polynomes : " << Optimal_Nb_Simulation_vdc_10 << std::endl;
//
//	delete vdc_gen_10;
//	delete vdc_gen_alias_10;
//	delete ngen_vdc_10;
//	delete vectorG_vdc_10;
//	delete BSamerican_vdc_10;
//
//	////////////////// 1.1.3) Vanille (Sobol - Laguerre Polynomes n = 6) //////////////////////////////////////
//
//	UniformGenerator* sobol_ugen_10 = new Sobol();
//	Normal* ngen_sobol_10 = new NormalBoxMuller(sobol_ugen_10, 0., 1.);
//
//	GaussianVector* vectorG_sobol_10 = new GaussianVectorCholesky(ngen_sobol_10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_10 = new BSEulerND(vectorG_sobol_10, spot_m, rate);
//
//	AmericanMonteCarlo_basket USMC_sobol_10(N, bsktcall, BSamerican_sobol_10, basefunc_Laguerre_6);
//
//	clock::time_point start_US_sobol_10 = clock::now();
//	USMC_sobol_10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_10 = clock::now();
//	clock::duration execution_timeUS_sobol_10 = end_US_sobol_10 - start_US_sobol_10;
//	double price_US_sobol_10 = USMC_sobol_10.GetPrice();
//	double var_US_sobol_10 = USMC_sobol_10.GetVariance();
//
//	USMC_sobol_10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_10 = USMC_sobol_10.GetNbSimul();
//
//	std::cout << " US Case - exec time for vanilla MonteCarlo with Sobol sequence and 6 order of Laguerre polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_10).count() << std::endl;
//	std::cout << " US Case - price for vanilla MonteCarlo with Sobol sequence and 6 order of Laguerre polynomes : " << price_US_sobol_10 << std::endl;
//	std::cout << " US Case - variance for vanilla MonteCarlo  with Sobol sequence and 6 order of Laguerre polynomes : " << var_US_sobol_10 << std::endl;
//	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo  with Sobol sequence : and 6 order of Laguerre polynomes " << Optimal_Nb_Simulation_sobol_10 << std::endl;
//
//	delete sobol_ugen_10;
//	delete ngen_sobol_10;
//	delete vectorG_sobol_10;
//	delete BSamerican_sobol_10; 
//
//
/////////////////////////////////////// VARIANCE REDUCTION TECHNIQUES ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//	//////////////////////////// CONTROL VARIABLE /////////////////////////////////////////////////////////////////////
//
//
//
//	///create both the control variate payoff and its close_formula price in order to proceed to the contrl variate technique
//
//	PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(W, spot_m, 100.);
//
//	ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(W, spot_m, CovarMatrix, K,
//		rate, endTime);
//
//	double df = exp(-rate * endTime);
//	double price_cf = CFbaskt->operator()(spot_m, df);
//
//
//	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// 2.4) Controle variable (Ecuyer Combined - Laguerre Polynomes n = 6) //////////////////////////////////////
//
//	UniformGenerator* ugen_VC10 = new EcuyerCombined();
//	Normal* ngen_VC10 = new NormalBoxMuller(ugen_VC10, 0., 1.);
//
//	GaussianVector* vectorG_VC10 = new GaussianVectorCholesky(ngen_VC10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_VC10 = new BSEulerND(vectorG_VC10, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_VC10(N, bsktcall, bsktcallCV, BSamerican_VC10,
//		basefunc_Laguerre_6, price_cf);
//
//	clock::time_point start_US_VC10 = clock::now();
//	USMC_VC10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_VC10 = clock::now();
//	clock::duration execution_timeUS_VC10 = end_US_VC10 - start_US_VC10;
//	double price_US_VC10 = USMC_VC10.GetPrice();
//	double var_US_VC10 = USMC_VC10.GetVariance();
//
//	USMC_VC10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_VC10 = USMC_VC10.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo and 6 order of Laguerre polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_VC10).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo and and 6 order of Laguerre polynomes : " << price_US_VC10 << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo and 6 order of Laguerre polynomes : " << var_US_VC10 << std::endl;
//	std::cout << " US Case - optimal number of simulations for control variable MonteCarlo and 6 order of Laguerre polynomes : " << Optimal_Nb_Simulation_VC10 << std::endl;
//
//	delete ugen_VC10;
//	delete ngen_VC10;
//	delete vectorG_VC10;
//	delete BSamerican_VC10;
//
//	////////////////// 2.4.2) Control variable (VDC - Laguerre Polynomes n = 6) //////////////////////////////////////
//
//	UniformGenerator* vdc_gen_VC10 = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_VC10 = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_VC10 = new NormalBoxMullerVDC(vdc_gen_VC10, vdc_gen_alias_VC10, 0., 1.);
//
//	GaussianVector* vectorG_vdc_VC10 = new GaussianVectorCholesky(ngen_vdc_VC10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_VC10 = new BSEulerND(vectorG_vdc_VC10, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_VC_vdc_Laguerre10(N, bsktcall, bsktcallCV, BSamerican_vdc_VC10,
//		basefunc_Laguerre_6, price_cf);
//
//	clock::time_point start_US_vdc_VC10 = clock::now();
//	USMC_VC_vdc_Laguerre10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_VC10 = clock::now();
//	clock::duration execution_timeUS_vdc_VC10 = end_US_vdc_VC10 - start_US_vdc_VC10;
//	double price_US_vdc_VC10 = USMC_VC_vdc_Laguerre10.GetPrice();
//	double var_US_vdc_VC10 = USMC_VC_vdc_Laguerre10.GetVariance();
//
//	USMC_VC_vdc_Laguerre10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_VC10 = USMC_VC_vdc_Laguerre10.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_VC10).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynomes : " << price_US_vdc_VC10 << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo  with VDC sequence and 6 order of Laguerre polynomes : " << var_US_vdc_VC10 << std::endl;
//	std::cout << " US Case - optimal number of simulations for control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynomes : " << Optimal_Nb_Simulation_vdc_VC10 << std::endl;
//
//	delete vdc_gen_VC10;
//	delete vdc_gen_alias_VC10;
//	delete ngen_vdc_VC10;
//	delete vectorG_vdc_VC10;
//	delete BSamerican_vdc_VC10;
//
//	////////////////// 2.4.3) Control variable (Sobol - Laguerre Polynomes n = 6) //////////////////////////////////////
//
//	UniformGenerator* sobol_ugen_CV10 = new Sobol();
//	Normal* ngen_sobol_CV10 = new NormalBoxMuller(sobol_ugen_CV10, 0., 1.);
//
//	GaussianVector* vectorG_sobol_CV10 = new GaussianVectorCholesky(ngen_sobol_CV10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_CV10 = new BSEulerND(vectorG_sobol_CV10, spot_m, rate);
//
//	AmericanMonteCarlo_basket_controlevariable USMC_sobol_CV_Laguerre10(N, bsktcall, bsktcallCV, BSamerican_sobol_CV10, basefunc_Laguerre_6,price_cf);
//
//	clock::time_point start_US_sobol_CV10 = clock::now();
//	USMC_sobol_CV_Laguerre10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_CV10 = clock::now();
//	clock::duration execution_timeUS_sobol_CV10 = end_US_sobol_CV10 - start_US_sobol_CV10;
//	double price_US_sobol_CV10 = USMC_sobol_CV_Laguerre10.GetPrice();
//	double var_US_sobol_CV10 = USMC_sobol_CV_Laguerre10.GetVariance();
//
//	USMC_sobol_CV_Laguerre10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_CV10 = USMC_sobol_CV_Laguerre10.GetNbSimul();
//
//	std::cout << " US Case - exec time for control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynomes  : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_CV10).count() << std::endl;
//	std::cout << " US Case - price for control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynomes : " << price_US_sobol_CV10 << std::endl;
//	std::cout << " US Case - variance for control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynomes : " << var_US_sobol_CV10 << std::endl;
//	std::cout << " US Case - optimal number of simulations forcontrol variable  MonteCarlo  with Sobol sequence and 6 order of Laguerre polynomes : " << Optimal_Nb_Simulation_sobol_CV10 << std::endl;
//
//	delete sobol_ugen_CV10;
//	delete ngen_sobol_CV10;
//	delete vectorG_sobol_CV10;
//	delete BSamerican_sobol_CV10;
//
//
//	//////////////////////// 3.4) Antithetic (Ecuyer Combined - Laguerre polynome n = 6) ////////////////////////////
//
//	UniformGenerator* ugen_anti10 = new EcuyerCombined();
//	Normal* ngen_anti10 = new NormalBoxMuller(ugen_anti10, 0., 1.);
//
//	GaussianVector* vectorG_anti10 = new GaussianVectorCholesky(ngen_anti10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_anti10 = new BSEulerNDAntithetic(vectorG_anti10, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_anti_Laguerre10(N, bsktcall, BSamerican_anti10,
//		basefunc_Laguerre_6);
//
//	clock::time_point start_US_anti10 = clock::now();
//	USMC_anti_Laguerre10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_anti10 = clock::now();
//	clock::duration execution_timeUS_anti10 = end_US_anti10 - start_US_anti10;
//	double price_US_anti10 = USMC_anti_Laguerre10.GetPrice();
//	double var_US_anti10 = USMC_anti_Laguerre10.GetVariance();
//
//	USMC_anti_Laguerre10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_anti10 = USMC_anti_Laguerre10.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with 6 order of Laguerre polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti10).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with  6 order of Laguerre polynome : " << price_US_anti10 << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with  6 order of Laguerre polynome: " << var_US_anti10 << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with  6 order of Laguerre polynome : " << Optimal_Nb_Simulation_anti10 << std::endl;
//
//	delete ugen_anti10;
//	delete ngen_anti10;
//	delete vectorG_anti10;
//	delete BSamerican_anti10;
//
//
//	//////////////////////// 3.4.2) Antithetic (VDC - Laguerre polynome n = 6) ////////////////////////////
//
//	UniformGenerator* vdc_gen_anti10 = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_anti10 = new VanDerCorput(2, 1);
//	Normal* ngen_vdc_anti10 = new NormalBoxMullerVDC(vdc_gen_anti10, vdc_gen_alias_anti10, 0., 1.);
//
//	GaussianVector* vectorG_vdc_anti10 = new GaussianVectorCholesky(ngen_vdc_anti10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_anti10 = new BSEulerNDAntithetic(vectorG_vdc_anti10, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_vdc_anti_Laguerre10(N, bsktcall, BSamerican_vdc_anti10,
//		basefunc_Laguerre_6);
//
//	clock::time_point start_US_vdc_anti10 = clock::now();
//	USMC_vdc_anti_Laguerre10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_anti10 = clock::now();
//	clock::duration execution_timeUS_vdc_anti10 = end_US_vdc_anti10 - start_US_vdc_anti10;
//	double price_US_vdc_anti10 = USMC_vdc_anti_Laguerre10.GetPrice();
//	double var_US_vdc_anti10 = USMC_vdc_anti_Laguerre10.GetVariance();
//
//	USMC_vdc_anti_Laguerre10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_vdc_anti10 = USMC_vdc_anti_Laguerre10.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with VDC sequence and with  6 order of Laguerre polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti10).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with VDC sequence and with  6 order of Laguerre polynome : " << price_US_vdc_anti10 << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with VDC sequence and with  6 order of Laguerre polynome : " << var_US_vdc_anti10 << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with VDC sequence and with  6 order of Laguerre polynome : " << Optimal_Nb_Simulation_vdc_anti10 << std::endl;
//
//	delete vdc_gen_anti10;
//	delete vdc_gen_alias_anti10;
//	delete ngen_vdc_anti10;
//	delete vectorG_vdc_anti10;
//	delete BSamerican_vdc_anti10;
//
//	//////////////////////// 3.4.3) Antithetic (Sobol - Laguerre Polynome n = 6) ////////////////////////////
//
//	UniformGenerator* gen_sobol_anti10 = new Sobol();
//	Normal* ngen_sobol_anti10 = new NormalBoxMuller(gen_sobol_anti10, 0., 1.);
//
//	GaussianVector* vectorG_sobol_anti10 = new GaussianVectorCholesky(ngen_sobol_anti10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_anti10 = new BSEulerNDAntithetic(vectorG_sobol_anti10, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic USMC_sobol_anti_Laguerre10(N, bsktcall, BSamerican_sobol_anti10,
//		basefunc_Laguerre_6);
//
//	clock::time_point start_US_sobol_anti10 = clock::now();
//	USMC_sobol_anti_Laguerre10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_anti10 = clock::now();
//	clock::duration execution_timeUS_sobol_anti10 = end_US_sobol_anti10 - start_US_sobol_anti10;
//	double price_US_sobol_anti10 = USMC_sobol_anti_Laguerre10.GetPrice();
//	double var_US_sobol_anti10 = USMC_sobol_anti_Laguerre10.GetVariance();
//
//	USMC_sobol_anti_Laguerre10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_sobol_anti10 = USMC_sobol_anti_Laguerre10.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic MonteCarlo with Sobol sequence and  with  6 order of Laguerre polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti10).count() << std::endl;
//	std::cout << " US Case - price for antithetic MonteCarlo with Sobol sequence and with  6 order of Laguerre polynome : " << price_US_sobol_anti10 << std::endl;
//	std::cout << " US Case - variance for antithetic MonteCarlo with Sobol sequence and with  6 order of Laguerre polynome : " << var_US_sobol_anti10 << std::endl;
//	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with Sobol sequence and with  6 order of Laguerre polynome : " << Optimal_Nb_Simulation_sobol_anti10 << std::endl;
//
//	delete gen_sobol_anti10;
//	delete ngen_sobol_anti10;
//	delete vectorG_sobol_anti10;
//	delete BSamerican_sobol_anti10;
//
//	//////////////////////// 4.4) VC + Antithetic (Ecuyer Combined - Laguerre polynome n = 10) ////////////////////////////
//
//	UniformGenerator* ugen_anti_cv10 = new EcuyerCombined();
//	Normal* ngen_anti_cv10 = new NormalBoxMuller(ugen_anti_cv10, 0., 1.);
//
//	GaussianVector* vectorG_anti_cv10 = new GaussianVectorCholesky(ngen_anti_cv10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_anti_cv10 = new BSEulerNDAntithetic(vectorG_anti_cv10, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_anti_CV_Laguerre10(N, bsktcall, bsktcallCV, BSamerican_anti_cv10, basefunc_Laguerre_6, price_cf);
//
//
//	clock::time_point start_US_anti_cv10 = clock::now();
//	USMC_anti_CV_Laguerre10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_anti_cv10 = clock::now();
//	clock::duration execution_timeUS_anti_cv10 = end_US_anti_cv10 - start_US_anti_cv10;
//	double price_US_anti_cv10 = USMC_anti_CV_Laguerre10.GetPrice();
//	double var_US_anti_cv10 = USMC_anti_CV_Laguerre10.GetVariance();
//
//	USMC_anti_CV_Laguerre10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_Nb_Simulation_anti_cv10 = USMC_anti_CV_Laguerre10.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo and 6 order of Laguerre polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti_cv10).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo and 6 order of Laguerre polynome : " << price_US_anti_cv10 << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo and 6 order of Laguerre polynome : " << var_US_anti_cv10 << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo and 6 order of Laguerre polynome : " << Optimal_Nb_Simulation_anti_cv10 << std::endl;
//
//	delete ugen_anti_cv10;
//	delete ngen_anti_cv10;
//	delete vectorG_anti_cv10;
//	delete BSamerican_anti_cv10;
//
//	//////////////////////// 4.4.2) VC + Antithetic (VDC - Laguerre polynome n = 10) ////////////////////////////
//
//	UniformGenerator* vdc_gen_anti_cv10 = new VanDerCorput(2, 1);
//	UniformGenerator* vdc_gen_alias_anti_cv10 = new VanDerCorput(11, 1);
//	Normal* ngen_vdc_anti_cv10 = new NormalBoxMullerVDC(vdc_gen_anti_cv10, vdc_gen_alias_anti_cv10, 0., 1.);
//
//	GaussianVector* vectorG_vdc_anti_cv10 = new GaussianVectorCholesky(ngen_vdc_anti_cv10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_vdc_anti_cv10 = new BSEulerNDAntithetic(vectorG_vdc_anti_cv10, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_vdc_anti_CV_Laguerre10(N, bsktcall, bsktcallCV, BSamerican_vdc_anti_cv10, basefunc_Laguerre_6, price_cf);
//
//
//	clock::time_point start_US_vdc_anti_cv10 = clock::now();
//	USMC_vdc_anti_CV_Laguerre10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_vdc_anti_cv10 = clock::now();
//	clock::duration execution_timeUS_vdc_anti_cv10 = end_US_vdc_anti_cv10 - start_US_vdc_anti_cv10;
//	double price_US_vdc_anti_cv10 = USMC_vdc_anti_CV_Laguerre10.GetPrice();
//	double var_US_vdc_anti_cv10 = USMC_vdc_anti_CV_Laguerre10.GetVariance();
//
//	USMC_vdc_anti_CV_Laguerre10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_vdc_Nb_Simulation_anti_cv10 = USMC_vdc_anti_CV_Laguerre10.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti_cv10).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo  with VDC sequence and 6 order of Laguerre polynome : " << price_US_vdc_anti_cv10 << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynome : " << var_US_vdc_anti_cv10 << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynome : " << Optimal_vdc_Nb_Simulation_anti_cv10 << std::endl;
//
//	delete vdc_gen_anti_cv10;
//	delete vdc_gen_alias_anti_cv10;
//	delete ngen_vdc_anti_cv10;
//	delete vectorG_vdc_anti_cv10;
//	delete BSamerican_vdc_anti_cv10;
//
//	//////////////////////// 4.4.3) VC + Antithetic (Sobol - Laguerre polynome n = 10) ////////////////////////////
//
//	UniformGenerator* sobol_gen_anti_cv10 = new Sobol();
//	Normal* ngen_sobol_anti_cv10 = new NormalBoxMuller(sobol_gen_anti_cv10, 0., 1.);
//
//	GaussianVector* vectorG_sobol_anti_cv10 = new GaussianVectorCholesky(ngen_sobol_anti_cv10, Sigma, Correl, CovarMatrix);
//	RandomProcess* BSamerican_sobol_anti_cv10 = new BSEulerNDAntithetic(vectorG_sobol_anti_cv10, spot_m, rate);
//
//	AmericanMonteCarlo_basket_Antithetic_CV USMC_sobol_anti_CV_Laguerre10(N, bsktcall, bsktcallCV, BSamerican_sobol_anti_cv10, basefunc_Laguerre_6, price_cf);
//
//
//	clock::time_point start_US_sobol_anti_cv10 = clock::now();
//	USMC_sobol_anti_CV_Laguerre10.Simulate(startTime, endTime, nbsteps);
//	clock::time_point end_US_sobol_anti_cv10 = clock::now();
//	clock::duration execution_timeUS_sobol_anti_cv10 = end_US_sobol_anti_cv10 - start_US_sobol_anti_cv10;
//	double price_US_sobol_anti_cv10 = USMC_sobol_anti_CV_Laguerre10.GetPrice();
//	double var_US_sobol_anti_cv10 = USMC_sobol_anti_CV_Laguerre10.GetVariance();
//
//	USMC_sobol_anti_CV_Laguerre10.OptimalNbSimul(tolerated_error);
//	size_t Optimal_sobol_Nb_Simulation_anti_cv10 = USMC_sobol_anti_CV_Laguerre10.GetNbSimul();
//
//	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynome: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti_cv10).count() << std::endl;
//	std::cout << " US Case - price for antithetic and control variable MonteCarlo  with Sobol sequence and 6 order of Laguerre polynome : " << price_US_sobol_anti_cv10 << std::endl;
//	std::cout << " US Case - variance for antithetic and control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynome : " << var_US_sobol_anti_cv10 << std::endl;
//	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynome : " << Optimal_sobol_Nb_Simulation_anti_cv10 << std::endl;
//
//	delete sobol_gen_anti_cv10;
//	delete ngen_sobol_anti_cv10;
//	delete vectorG_sobol_anti_cv10;
//	delete BSamerican_sobol_anti_cv10;

	return 0;
}