#include "sde.hpp"
#include "MonteCarlo.h"
#include <iostream>
#include <fstream>
#include <chrono> 

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

/*
	std::cout<<"LowerDiag Cholesky"<<std::endl;
	matrix LowerDiag = CovarMatrix.Cholesky();
	LowerDiag.Print();
	
	matrix v(CovarMatrix.nb_rows(),CovarMatrix.nb_rows());
	matrix d(CovarMatrix.nb_rows(),1);
	int n = CovarMatrix.nb_rows();
	int nrot= 0;
	jacobi(CovarMatrix,n,d, v, nrot);
	
	std::cout<<"CovarMatrix after function called"<<std::endl;
	CovarMatrix.Print();
	std::cout<<"d matrix of eigevalues"<<std::endl;
	d.Print();
	std::cout<<"v matrix of eigenvectors"<<std::endl;
	v.Print();
	std::cout<<"nrot"<<nrot<<std::endl;
	
	d.Diagonalization();
	std::cout<<"d matrix of eigevalues"<<std::endl;
	d.Print();
	std::cout<<"d sqrt"<<std::endl;
	d.SQRT();
	d.Print();
	matrix l = v*d;
	std::cout<<"output matrix"<<std::endl;
	l.Print();
*/

/*
	UniformGenerator* ugen = new EcuyerCombined();
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);

	GaussianVector* corrG = new GaussianVectorDiag(ngen, Sigma, Correl, CovarMatrix);
	GaussianVector* corrG2 = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	matrix t = corrG->GetB();
	matrix t2 = corrG2->GetB();
	std::cout<<"output matrix"<<std::endl;
	t.Print();
	t2.Print();
	
	delete ugen;
	delete ngen;
	delete corrG;
	delete corrG2;
*/
/////////////////////////////////////////////////////////////////////////////////////

//1. NO VARIANCE REDUCTION

	//RANDOM NUMBER GENERATION
	UniformGenerator* ugen = new EcuyerCombined();
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	//RANDOM PROCESS
	RandomProcess* path = new BSEulerND(corrG,spot_m,rate);

	//PAYOFF
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


/////////////////////////////////////////////////////////////////////////////////////
/*
//2.QUASI MC WITH SOBOL SEQUENCE
	//RANDOM NUMBER GENERATION
	UniformGenerator* ugen = new Sobol();
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	//RANDOM PROCESS
	RandomProcess* path = new BSEulerND(corrG,spot_m,rate);

	//PAYOFF
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
	std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;

*/
/////////////////////////////////////////////////////////////////////////////////////

//3.ANTHITETIC MC
	//RANDOM NUMBER GENERATION
	UniformGenerator* ugen = new EcuyerCombined();
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	//RANDOM PROCESS
	RandomProcess* path = new BSEulerNDAntithetic(corrG,spot_m,rate);

	//PAYOFF
	PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);

	//MONTE CARLO
	double tolerated_error = 0.01;
	size_t CompMC_Nb_Simulation = 2000;
	
	EuropeanBasket_Antithetic MC(CompMC_Nb_Simulation, bsktcall, path);
	
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
	std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;

/////////////////////////////////////////////////////////////////////////////////////

//4.CV MC
	//RANDOM NUMBER GENERATION
	
	UniformGenerator* ugen = new EcuyerCombined();
	
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	// UniformGenerator* ugen = new Sobol();
	
	//RANDOM PROCESS
	RandomProcess* path = new BSEulerND(corrG,spot_m,rate);

	//PAYOFF
	PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);
	
	PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(Weights, spot_m,K);
	ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall
	(Weights, spot_m, CovarMatrix, K,rate, endTime);
	double CV_closedprice = CFbaskt->operator()(spot_m,df);
	
	//MONTE CARLO
	double tolerated_error = 0.01;
	size_t CompMC_Nb_Simulation = 2000;
	
	EuropeanBasket_controlvariable MC(CompMC_Nb_Simulation,
	bsktcall,bsktcallCV, path,CV_closedprice);
	
	//First Run: Companion MC to determinate the optimal number of simulation
	
	clock::time_point startComp = clock::now();
	MC.Simulate(startTime,endTime,nbsteps);
	clock::time_point endComp = clock::now();
	clock::duration execution_timeComp = endComp - startComp;
	
	//dynamics.Simulate(startTime, endTime, nbsteps);
	//matrix chemin = dynamics.GetAllPaths();
	//chemin.Print();

	UniformGenerator* ugen2 = new EcuyerCombined();
	RandomGenerator* ngen2 = new NormalBoxMuller(ugen2, 0., 1.);
	// RandomProcess* chemin = new BSEuler1D(ngen2, spot, rate, vol);
/*
	
	MC.OptimalNbSimul(tolerated_error);
	
	size_t Optimal_Nb_Simulation = MC.GetNbSimul();
	
	std::cout << "NB SIMULATION COMP: " << CompMC_Nb_Simulation << std::endl;
	std::cout << "ERROR TOLERATED: " << tolerated_error << std::endl;
	std::cout << "PRIX COMP: " << MC.GetPrice() << std::endl;
	std::cout << "Variance COMP: " << MC.GetVariance() << std::endl;
	std::cout << "NB SIMULATION OPTI: " << Optimal_Nb_Simulation << std::endl;
	std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;
	delete CFbaskt;
	delete bsktcallCV;


/////////////////////////////////////////////////////////////////////////////////////

/*
//5.ANTHITETIC + SOBOL MC
	//RANDOM NUMBER GENERATION
	UniformGenerator* ugen = new Sobol();
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	//RANDOM PROCESS
	RandomProcess* path = new BSEulerNDAntithetic(corrG,spot_m,rate);

	//PAYOFF
	PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);

	//MONTE CARLO
	double tolerated_error = 0.01;
	size_t CompMC_Nb_Simulation = 2000;
	
	EuropeanBasket_Antithetic MC(CompMC_Nb_Simulation, bsktcall, path);
	
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
	std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;
*/
/////////////////////////////////////////////////////////////////////////////////////
/*
//6. QUASI + CV MC
	//RANDOM NUMBER GENERATION
	UniformGenerator* ugen = new Sobol();
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	GaussianVector* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	//RANDOM PROCESS
	RandomProcess* path = new BSEulerND(corrG,spot_m,rate);

	//PAYOFF
	PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);
	PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(Weights, spot_m,K);
	ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(Weights, spot_m, CovarMatrix, K,rate, endTime);

	double CV_closedprice = CFbaskt->operator()(spot_m,df);
	// std::cout << "closed formula for the basket option" << CV_closedprice << std::endl;
	//MONTE CARLO
	double tolerated_error = 0.01;
	size_t CompMC_Nb_Simulation = 2000;
	
	EuropeanBasket_controlvariable MC(CompMC_Nb_Simulation, bsktcall,bsktcallCV, path,CV_closedprice);
	
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
	std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;
	delete CFbaskt;
	delete bsktcallCV;

*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//7. Anti + CV MC
	//RANDOM NUMBER GENERATION
	UniformGenerator* ugen = new EcuyerCombined();
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	GaussianVectorCholesky* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	//RANDOM PROCESS
	BSEulerNDAntithetic* path = new BSEulerNDAntithetic(corrG,spot_m,rate);

	//PAYOFF
	PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);
	PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(Weights, spot_m,K);
	ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(Weights, spot_m, CovarMatrix, K,rate, endTime);

	double CV_closedprice = CFbaskt->operator()(spot_m,df);
	// std::cout << "closed formula for the basket option" << CV_closedprice << std::endl;
	//MONTE CARLO
	double tolerated_error = 0.01;
	size_t CompMC_Nb_Simulation = 2000;
	
	EuropeanBasket_Antithetic_CV MC(CompMC_Nb_Simulation,
	bsktcall,bsktcallCV, path,CV_closedprice);
	
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
	std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;
	delete CFbaskt;
	delete bsktcallCV;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//8. 3 TECHNIQUES USED
	//RANDOM NUMBER GENERATION
	UniformGenerator* ugen = new Sobol();
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	GaussianVectorCholesky* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	
	//RANDOM PROCESS
	BSEulerNDAntithetic* path = new BSEulerNDAntithetic(corrG,spot_m,rate);

	//PAYOFF
	PayOffBasket* bsktcall = new PayOffBasketCall(Weights, spot_m,K);
	PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(Weights, spot_m,K);
	ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(Weights, spot_m, CovarMatrix, K,rate, endTime);

	double CV_closedprice = CFbaskt->operator()(spot_m,df);
	// std::cout << "closed formula for the basket option" << CV_closedprice << std::endl;
	//MONTE CARLO
	double tolerated_error = 0.01;
	size_t CompMC_Nb_Simulation = 2000;
	
	EuropeanBasket_Antithetic_CV MC(CompMC_Nb_Simulation, bsktcall,bsktcallCV, path,CV_closedprice);
	
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
	std::cout << "EXECUTION TIME COMP: " << std::chrono::duration <double,std::ratio<1>> (execution_timeComp).count() << std::endl;
	delete CFbaskt;
	delete bsktcallCV;

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

	delete ugen;
	delete ngen;
	delete corrG;
	delete path;
	delete bsktcall;

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

	return 0;
}
