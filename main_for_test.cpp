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
	double K = 100;
	double rate = 0.05;
	double dt = (endTime - startTime) / nbsteps;

	double tolerated_error = 0.01;

	size_t N = 2000;

	std::vector<std::vector<double>> Spot_vector = { {95.},{120.},{80.} };

	std::vector<std::vector<double>> Sigma_vector = { {0.25},{0.20},{0.15} };

	std::vector<std::vector<double>> Correl_mat = { {1.,-0.2,0.4},
													{-0.2,1.,0.4},
												   {0.4,0.4,1.} };

	std::vector<std::vector<double>> Weights_m = { {0.3,0.5,0.2} };

	matrix spot_m(Spot_vector);
	matrix Sigma(Sigma_vector);
	matrix Correl(Correl_mat);
	matrix W(Weights_m);
	matrix CovarMatrix = VarCovarMatrix(Sigma, Correl);


	//// need to choose the polynome basis with which one will approxime the continuation value of the option
	std::vector<basis_functions*> basefunc_Laguerre_6;

	basefunc_Laguerre_6.push_back(new Poly_Laguerre(0));
	basefunc_Laguerre_6.push_back(new Poly_Laguerre(1));
	basefunc_Laguerre_6.push_back(new Poly_Laguerre(2));
	basefunc_Laguerre_6.push_back(new Poly_Laguerre(3));
	basefunc_Laguerre_6.push_back(new Poly_Laguerre(4));
	basefunc_Laguerre_6.push_back(new Poly_Laguerre(5));
	basefunc_Laguerre_6.push_back(new Poly_Laguerre(6));

	//Create the payoff of the options 

	PayOffBasket* bsktcall = new PayOffBasketCall(W, spot_m, 100.);

	/////////////////// NO VARIANCE REDUCTION /////////////////////////////


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////// 1.4) Vanille (Ecuyer Combined - Laguerre Polynomes n = 6) //////////////////////////////////////

	UniformGenerator* ugen10 = new EcuyerCombined();
	Normal* ngen10 = new NormalBoxMuller(ugen10, 0., 1.);

	GaussianVector* vectorG10 = new GaussianVectorCholesky(ngen10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican10 = new BSEulerND(vectorG10, spot_m, rate);

	AmericanMonteCarlo_basket USMC_10(N, bsktcall, BSamerican10, basefunc_Laguerre_6);

	clock::time_point start_US_10 = clock::now();
	USMC_10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_10 = clock::now();
	clock::duration execution_timeUS_10 = end_US_10 - start_US_10;
	double price_US_10 = USMC_10.GetPrice();
	double var_US_10 = USMC_10.GetVariance();

	USMC_10.OptimalNbSimul(tolerated_error);
	size_t Optimal_Nb_Simulation_10 = USMC_10.GetNbSimul();

	std::cout << " US Case - exec time for vanilla MonteCarlo with 6 order of Laguerre polynomes: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_10).count() << std::endl;
	std::cout << " US Case - price for vanilla MonteCarlo with 6 order of Laguerre polynomes : " << price_US_10 << std::endl;
	std::cout << " US Case - variance for vanilla MonteCarlo with 6 order of Laguerre polynomes : " << var_US_10 << std::endl;
	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo with 6 order of Laguerre polynomes : " << Optimal_Nb_Simulation_10 << std::endl;

	delete ugen10;
	delete ngen10;
	delete vectorG10;
	delete BSamerican10;

	////////////////// 1.4.2) Vanille (VDC - Laguerre Polynomes n = 6) //////////////////////////////////////

	UniformGenerator* vdc_gen_10 = new VanDerCorput(2, 1);
	UniformGenerator* vdc_gen_alias_10 = new VanDerCorput(11, 1);
	Normal* ngen_vdc_10 = new NormalBoxMullerVDC(vdc_gen_10, vdc_gen_alias_10, 0., 1.);

	GaussianVector* vectorG_vdc_10 = new GaussianVectorCholesky(ngen_vdc_10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_vdc_10 = new BSEulerND(vectorG_vdc_10, spot_m, rate);

	AmericanMonteCarlo_basket USMC_vdc_10(N, bsktcall, BSamerican_vdc_10, basefunc_Laguerre_6);

	clock::time_point start_US_vdc_10 = clock::now();
	USMC_vdc_10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_vdc_10 = clock::now();
	clock::duration execution_timeUS_vdc_10 = end_US_vdc_10 - start_US_vdc_10;
	double price_US_vdc_10 = USMC_vdc_10.GetPrice();
	double var_US_vdc_10 = USMC_vdc_10.GetVariance();

	USMC_vdc_10.OptimalNbSimul(tolerated_error);
	size_t Optimal_Nb_Simulation_vdc_10 = USMC_vdc_10.GetNbSimul();

	std::cout << " US Case - exec time for vanilla MonteCarlo with VDC sequence and 6 order of Laguerre polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_10).count() << std::endl;
	std::cout << " US Case - price for vanilla MonteCarlo with VDC sequence and 6 order of Laguerre polynomes : " << price_US_vdc_10 << std::endl;
	std::cout << " US Case - variance for vanilla MonteCarlo  with VDC sequence and 6 order of Laguerre polynomes : " << var_US_vdc_10 << std::endl;
	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo  with VDC sequence and 6 order of Laguerre polynomes : " << Optimal_Nb_Simulation_vdc_10 << std::endl;

	delete vdc_gen_10;
	delete vdc_gen_alias_10;
	delete ngen_vdc_10;
	delete vectorG_vdc_10;
	delete BSamerican_vdc_10;

	////////////////// 1.1.3) Vanille (Sobol - Laguerre Polynomes n = 6) //////////////////////////////////////

	UniformGenerator* sobol_ugen_10 = new Sobol();
	Normal* ngen_sobol_10 = new NormalBoxMuller(sobol_ugen_10, 0., 1.);

	GaussianVector* vectorG_sobol_10 = new GaussianVectorCholesky(ngen_sobol_10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_sobol_10 = new BSEulerND(vectorG_sobol_10, spot_m, rate);

	AmericanMonteCarlo_basket USMC_sobol_10(N, bsktcall, BSamerican_sobol_10, basefunc_Laguerre_6);

	clock::time_point start_US_sobol_10 = clock::now();
	USMC_sobol_10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_sobol_10 = clock::now();
	clock::duration execution_timeUS_sobol_10 = end_US_sobol_10 - start_US_sobol_10;
	double price_US_sobol_10 = USMC_sobol_10.GetPrice();
	double var_US_sobol_10 = USMC_sobol_10.GetVariance();

	USMC_sobol_10.OptimalNbSimul(tolerated_error);
	size_t Optimal_Nb_Simulation_sobol_10 = USMC_sobol_10.GetNbSimul();

	std::cout << " US Case - exec time for vanilla MonteCarlo with Sobol sequence and 6 order of Laguerre polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_10).count() << std::endl;
	std::cout << " US Case - price for vanilla MonteCarlo with Sobol sequence and 6 order of Laguerre polynomes : " << price_US_sobol_10 << std::endl;
	std::cout << " US Case - variance for vanilla MonteCarlo  with Sobol sequence and 6 order of Laguerre polynomes : " << var_US_sobol_10 << std::endl;
	std::cout << " US Case - optimal number of simulations for vanilla MonteCarlo  with Sobol sequence : and 6 order of Laguerre polynomes " << Optimal_Nb_Simulation_sobol_10 << std::endl;

	delete sobol_ugen_10;
	delete ngen_sobol_10;
	delete vectorG_sobol_10;
	delete BSamerican_sobol_10; 


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


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// 2.4) Controle variable (Ecuyer Combined - Laguerre Polynomes n = 6) //////////////////////////////////////

	UniformGenerator* ugen_VC10 = new EcuyerCombined();
	Normal* ngen_VC10 = new NormalBoxMuller(ugen_VC10, 0., 1.);

	GaussianVector* vectorG_VC10 = new GaussianVectorCholesky(ngen_VC10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_VC10 = new BSEulerND(vectorG_VC10, spot_m, rate);

	AmericanMonteCarlo_basket_controlevariable USMC_VC10(N, bsktcall, bsktcallCV, BSamerican_VC10,
		basefunc_Laguerre_6, price_cf);

	clock::time_point start_US_VC10 = clock::now();
	USMC_VC10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_VC10 = clock::now();
	clock::duration execution_timeUS_VC10 = end_US_VC10 - start_US_VC10;
	double price_US_VC10 = USMC_VC10.GetPrice();
	double var_US_VC10 = USMC_VC10.GetVariance();

	USMC_VC10.OptimalNbSimul(tolerated_error);
	size_t Optimal_Nb_Simulation_VC10 = USMC_VC10.GetNbSimul();

	std::cout << " US Case - exec time for control variable MonteCarlo and 6 order of Laguerre polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_VC10).count() << std::endl;
	std::cout << " US Case - price for control variable MonteCarlo and and 6 order of Laguerre polynomes : " << price_US_VC10 << std::endl;
	std::cout << " US Case - variance for control variable MonteCarlo and 6 order of Laguerre polynomes : " << var_US_VC10 << std::endl;
	std::cout << " US Case - optimal number of simulations for control variable MonteCarlo and 6 order of Laguerre polynomes : " << Optimal_Nb_Simulation_VC10 << std::endl;

	delete ugen_VC10;
	delete ngen_VC10;
	delete vectorG_VC10;
	delete BSamerican_VC10;

	////////////////// 2.4.2) Control variable (VDC - Laguerre Polynomes n = 6) //////////////////////////////////////

	UniformGenerator* vdc_gen_VC10 = new VanDerCorput(2, 1);
	UniformGenerator* vdc_gen_alias_VC10 = new VanDerCorput(11, 1);
	Normal* ngen_vdc_VC10 = new NormalBoxMullerVDC(vdc_gen_VC10, vdc_gen_alias_VC10, 0., 1.);

	GaussianVector* vectorG_vdc_VC10 = new GaussianVectorCholesky(ngen_vdc_VC10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_vdc_VC10 = new BSEulerND(vectorG_vdc_VC10, spot_m, rate);

	AmericanMonteCarlo_basket_controlevariable USMC_VC_vdc_Laguerre10(N, bsktcall, bsktcallCV, BSamerican_vdc_VC10,
		basefunc_Laguerre_6, price_cf);

	clock::time_point start_US_vdc_VC10 = clock::now();
	USMC_VC_vdc_Laguerre10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_vdc_VC10 = clock::now();
	clock::duration execution_timeUS_vdc_VC10 = end_US_vdc_VC10 - start_US_vdc_VC10;
	double price_US_vdc_VC10 = USMC_VC_vdc_Laguerre10.GetPrice();
	double var_US_vdc_VC10 = USMC_VC_vdc_Laguerre10.GetVariance();

	USMC_VC_vdc_Laguerre10.OptimalNbSimul(tolerated_error);
	size_t Optimal_Nb_Simulation_vdc_VC10 = USMC_VC_vdc_Laguerre10.GetNbSimul();

	std::cout << " US Case - exec time for control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynomes : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_VC10).count() << std::endl;
	std::cout << " US Case - price for control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynomes : " << price_US_vdc_VC10 << std::endl;
	std::cout << " US Case - variance for control variable MonteCarlo  with VDC sequence and 6 order of Laguerre polynomes : " << var_US_vdc_VC10 << std::endl;
	std::cout << " US Case - optimal number of simulations for control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynomes : " << Optimal_Nb_Simulation_vdc_VC10 << std::endl;

	delete vdc_gen_VC10;
	delete vdc_gen_alias_VC10;
	delete ngen_vdc_VC10;
	delete vectorG_vdc_VC10;
	delete BSamerican_vdc_VC10;

	////////////////// 2.4.3) Control variable (Sobol - Laguerre Polynomes n = 6) //////////////////////////////////////

	UniformGenerator* sobol_ugen_CV10 = new Sobol();
	Normal* ngen_sobol_CV10 = new NormalBoxMuller(sobol_ugen_CV10, 0., 1.);

	GaussianVector* vectorG_sobol_CV10 = new GaussianVectorCholesky(ngen_sobol_CV10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_sobol_CV10 = new BSEulerND(vectorG_sobol_CV10, spot_m, rate);

	AmericanMonteCarlo_basket_controlevariable USMC_sobol_CV_Laguerre10(N, bsktcall, bsktcallCV, BSamerican_sobol_CV10, basefunc_Laguerre_6,price_cf);

	clock::time_point start_US_sobol_CV10 = clock::now();
	USMC_sobol_CV_Laguerre10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_sobol_CV10 = clock::now();
	clock::duration execution_timeUS_sobol_CV10 = end_US_sobol_CV10 - start_US_sobol_CV10;
	double price_US_sobol_CV10 = USMC_sobol_CV_Laguerre10.GetPrice();
	double var_US_sobol_CV10 = USMC_sobol_CV_Laguerre10.GetVariance();

	USMC_sobol_CV_Laguerre10.OptimalNbSimul(tolerated_error);
	size_t Optimal_Nb_Simulation_sobol_CV10 = USMC_sobol_CV_Laguerre10.GetNbSimul();

	std::cout << " US Case - exec time for control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynomes  : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_CV10).count() << std::endl;
	std::cout << " US Case - price for control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynomes : " << price_US_sobol_CV10 << std::endl;
	std::cout << " US Case - variance for control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynomes : " << var_US_sobol_CV10 << std::endl;
	std::cout << " US Case - optimal number of simulations forcontrol variable  MonteCarlo  with Sobol sequence and 6 order of Laguerre polynomes : " << Optimal_Nb_Simulation_sobol_CV10 << std::endl;

	delete sobol_ugen_CV10;
	delete ngen_sobol_CV10;
	delete vectorG_sobol_CV10;
	delete BSamerican_sobol_CV10;


	//////////////////////// 3.4) Antithetic (Ecuyer Combined - Laguerre polynome n = 6) ////////////////////////////

	UniformGenerator* ugen_anti10 = new EcuyerCombined();
	Normal* ngen_anti10 = new NormalBoxMuller(ugen_anti10, 0., 1.);

	GaussianVector* vectorG_anti10 = new GaussianVectorCholesky(ngen_anti10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_anti10 = new BSEulerNDAntithetic(vectorG_anti10, spot_m, rate);

	AmericanMonteCarlo_basket_Antithetic USMC_anti_Laguerre10(N, bsktcall, BSamerican_anti10,
		basefunc_Laguerre_6);

	clock::time_point start_US_anti10 = clock::now();
	USMC_anti_Laguerre10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_anti10 = clock::now();
	clock::duration execution_timeUS_anti10 = end_US_anti10 - start_US_anti10;
	double price_US_anti10 = USMC_anti_Laguerre10.GetPrice();
	double var_US_anti10 = USMC_anti_Laguerre10.GetVariance();

	USMC_anti_Laguerre10.OptimalNbSimul(tolerated_error);
	size_t Optimal_Nb_Simulation_anti10 = USMC_anti_Laguerre10.GetNbSimul();

	std::cout << " US Case - exec time for antithetic MonteCarlo with 6 order of Laguerre polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti10).count() << std::endl;
	std::cout << " US Case - price for antithetic MonteCarlo with  6 order of Laguerre polynome : " << price_US_anti10 << std::endl;
	std::cout << " US Case - variance for antithetic MonteCarlo with  6 order of Laguerre polynome: " << var_US_anti10 << std::endl;
	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with  6 order of Laguerre polynome : " << Optimal_Nb_Simulation_anti10 << std::endl;

	delete ugen_anti10;
	delete ngen_anti10;
	delete vectorG_anti10;
	delete BSamerican_anti10;


	//////////////////////// 3.4.2) Antithetic (VDC - Laguerre polynome n = 6) ////////////////////////////

	UniformGenerator* vdc_gen_anti10 = new VanDerCorput(2, 1);
	UniformGenerator* vdc_gen_alias_anti10 = new VanDerCorput(2, 1);
	Normal* ngen_vdc_anti10 = new NormalBoxMullerVDC(vdc_gen_anti10, vdc_gen_alias_anti10, 0., 1.);

	GaussianVector* vectorG_vdc_anti10 = new GaussianVectorCholesky(ngen_vdc_anti10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_vdc_anti10 = new BSEulerNDAntithetic(vectorG_vdc_anti10, spot_m, rate);

	AmericanMonteCarlo_basket_Antithetic USMC_vdc_anti_Laguerre10(N, bsktcall, BSamerican_vdc_anti10,
		basefunc_Laguerre_6);

	clock::time_point start_US_vdc_anti10 = clock::now();
	USMC_vdc_anti_Laguerre10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_vdc_anti10 = clock::now();
	clock::duration execution_timeUS_vdc_anti10 = end_US_vdc_anti10 - start_US_vdc_anti10;
	double price_US_vdc_anti10 = USMC_vdc_anti_Laguerre10.GetPrice();
	double var_US_vdc_anti10 = USMC_vdc_anti_Laguerre10.GetVariance();

	USMC_vdc_anti_Laguerre10.OptimalNbSimul(tolerated_error);
	size_t Optimal_Nb_Simulation_vdc_anti10 = USMC_vdc_anti_Laguerre10.GetNbSimul();

	std::cout << " US Case - exec time for antithetic MonteCarlo with VDC sequence and with  6 order of Laguerre polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti10).count() << std::endl;
	std::cout << " US Case - price for antithetic MonteCarlo with VDC sequence and with  6 order of Laguerre polynome : " << price_US_vdc_anti10 << std::endl;
	std::cout << " US Case - variance for antithetic MonteCarlo with VDC sequence and with  6 order of Laguerre polynome : " << var_US_vdc_anti10 << std::endl;
	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with VDC sequence and with  6 order of Laguerre polynome : " << Optimal_Nb_Simulation_vdc_anti10 << std::endl;

	delete vdc_gen_anti10;
	delete vdc_gen_alias_anti10;
	delete ngen_vdc_anti10;
	delete vectorG_vdc_anti10;
	delete BSamerican_vdc_anti10;

	//////////////////////// 3.4.3) Antithetic (Sobol - Laguerre Polynome n = 6) ////////////////////////////

	UniformGenerator* gen_sobol_anti10 = new Sobol();
	Normal* ngen_sobol_anti10 = new NormalBoxMuller(gen_sobol_anti10, 0., 1.);

	GaussianVector* vectorG_sobol_anti10 = new GaussianVectorCholesky(ngen_sobol_anti10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_sobol_anti10 = new BSEulerNDAntithetic(vectorG_sobol_anti10, spot_m, rate);

	AmericanMonteCarlo_basket_Antithetic USMC_sobol_anti_Laguerre10(N, bsktcall, BSamerican_sobol_anti10,
		basefunc_Laguerre_6);

	clock::time_point start_US_sobol_anti10 = clock::now();
	USMC_sobol_anti_Laguerre10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_sobol_anti10 = clock::now();
	clock::duration execution_timeUS_sobol_anti10 = end_US_sobol_anti10 - start_US_sobol_anti10;
	double price_US_sobol_anti10 = USMC_sobol_anti_Laguerre10.GetPrice();
	double var_US_sobol_anti10 = USMC_sobol_anti_Laguerre10.GetVariance();

	USMC_sobol_anti_Laguerre10.OptimalNbSimul(tolerated_error);
	size_t Optimal_Nb_Simulation_sobol_anti10 = USMC_sobol_anti_Laguerre10.GetNbSimul();

	std::cout << " US Case - exec time for antithetic MonteCarlo with Sobol sequence and  with  6 order of Laguerre polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti10).count() << std::endl;
	std::cout << " US Case - price for antithetic MonteCarlo with Sobol sequence and with  6 order of Laguerre polynome : " << price_US_sobol_anti10 << std::endl;
	std::cout << " US Case - variance for antithetic MonteCarlo with Sobol sequence and with  6 order of Laguerre polynome : " << var_US_sobol_anti10 << std::endl;
	std::cout << " US Case - optimal number of simulations for antithetic MonteCarlo with Sobol sequence and with  6 order of Laguerre polynome : " << Optimal_Nb_Simulation_sobol_anti10 << std::endl;

	delete gen_sobol_anti10;
	delete ngen_sobol_anti10;
	delete vectorG_sobol_anti10;
	delete BSamerican_sobol_anti10;

	//////////////////////// 4.4) VC + Antithetic (Ecuyer Combined - Laguerre polynome n = 10) ////////////////////////////

	UniformGenerator* ugen_anti_cv10 = new EcuyerCombined();
	Normal* ngen_anti_cv10 = new NormalBoxMuller(ugen_anti_cv10, 0., 1.);

	GaussianVector* vectorG_anti_cv10 = new GaussianVectorCholesky(ngen_anti_cv10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_anti_cv10 = new BSEulerNDAntithetic(vectorG_anti_cv10, spot_m, rate);

	AmericanMonteCarlo_basket_Antithetic_CV USMC_anti_CV_Laguerre10(N, bsktcall, bsktcallCV, BSamerican_anti_cv10, basefunc_Laguerre_6, price_cf);


	clock::time_point start_US_anti_cv10 = clock::now();
	USMC_anti_CV_Laguerre10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_anti_cv10 = clock::now();
	clock::duration execution_timeUS_anti_cv10 = end_US_anti_cv10 - start_US_anti_cv10;
	double price_US_anti_cv10 = USMC_anti_CV_Laguerre10.GetPrice();
	double var_US_anti_cv10 = USMC_anti_CV_Laguerre10.GetVariance();

	USMC_anti_CV_Laguerre10.OptimalNbSimul(tolerated_error);
	size_t Optimal_Nb_Simulation_anti_cv10 = USMC_anti_CV_Laguerre10.GetNbSimul();

	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo and 6 order of Laguerre polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti_cv10).count() << std::endl;
	std::cout << " US Case - price for antithetic and control variable MonteCarlo and 6 order of Laguerre polynome : " << price_US_anti_cv10 << std::endl;
	std::cout << " US Case - variance for antithetic and control variable MonteCarlo and 6 order of Laguerre polynome : " << var_US_anti_cv10 << std::endl;
	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo and 6 order of Laguerre polynome : " << Optimal_Nb_Simulation_anti_cv10 << std::endl;

	delete ugen_anti_cv10;
	delete ngen_anti_cv10;
	delete vectorG_anti_cv10;
	delete BSamerican_anti_cv10;

	//////////////////////// 4.4.2) VC + Antithetic (VDC - Laguerre polynome n = 10) ////////////////////////////

	UniformGenerator* vdc_gen_anti_cv10 = new VanDerCorput(2, 1);
	UniformGenerator* vdc_gen_alias_anti_cv10 = new VanDerCorput(11, 1);
	Normal* ngen_vdc_anti_cv10 = new NormalBoxMullerVDC(vdc_gen_anti_cv10, vdc_gen_alias_anti_cv10, 0., 1.);

	GaussianVector* vectorG_vdc_anti_cv10 = new GaussianVectorCholesky(ngen_vdc_anti_cv10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_vdc_anti_cv10 = new BSEulerNDAntithetic(vectorG_vdc_anti_cv10, spot_m, rate);

	AmericanMonteCarlo_basket_Antithetic_CV USMC_vdc_anti_CV_Laguerre10(N, bsktcall, bsktcallCV, BSamerican_vdc_anti_cv10, basefunc_Laguerre_6, price_cf);


	clock::time_point start_US_vdc_anti_cv10 = clock::now();
	USMC_vdc_anti_CV_Laguerre10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_vdc_anti_cv10 = clock::now();
	clock::duration execution_timeUS_vdc_anti_cv10 = end_US_vdc_anti_cv10 - start_US_vdc_anti_cv10;
	double price_US_vdc_anti_cv10 = USMC_vdc_anti_CV_Laguerre10.GetPrice();
	double var_US_vdc_anti_cv10 = USMC_vdc_anti_CV_Laguerre10.GetVariance();

	USMC_vdc_anti_CV_Laguerre10.OptimalNbSimul(tolerated_error);
	size_t Optimal_vdc_Nb_Simulation_anti_cv10 = USMC_vdc_anti_CV_Laguerre10.GetNbSimul();

	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynome : " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_vdc_anti_cv10).count() << std::endl;
	std::cout << " US Case - price for antithetic and control variable MonteCarlo  with VDC sequence and 6 order of Laguerre polynome : " << price_US_vdc_anti_cv10 << std::endl;
	std::cout << " US Case - variance for antithetic and control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynome : " << var_US_vdc_anti_cv10 << std::endl;
	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo with VDC sequence and 6 order of Laguerre polynome : " << Optimal_vdc_Nb_Simulation_anti_cv10 << std::endl;

	delete vdc_gen_anti_cv10;
	delete vdc_gen_alias_anti_cv10;
	delete ngen_vdc_anti_cv10;
	delete vectorG_vdc_anti_cv10;
	delete BSamerican_vdc_anti_cv10;

	//////////////////////// 4.4.3) VC + Antithetic (Sobol - Laguerre polynome n = 10) ////////////////////////////

	UniformGenerator* sobol_gen_anti_cv10 = new Sobol();
	Normal* ngen_sobol_anti_cv10 = new NormalBoxMuller(sobol_gen_anti_cv10, 0., 1.);

	GaussianVector* vectorG_sobol_anti_cv10 = new GaussianVectorCholesky(ngen_sobol_anti_cv10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_sobol_anti_cv10 = new BSEulerNDAntithetic(vectorG_sobol_anti_cv10, spot_m, rate);

	AmericanMonteCarlo_basket_Antithetic_CV USMC_sobol_anti_CV_Laguerre10(N, bsktcall, bsktcallCV, BSamerican_sobol_anti_cv10, basefunc_Laguerre_6, price_cf);


	clock::time_point start_US_sobol_anti_cv10 = clock::now();
	USMC_sobol_anti_CV_Laguerre10.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_sobol_anti_cv10 = clock::now();
	clock::duration execution_timeUS_sobol_anti_cv10 = end_US_sobol_anti_cv10 - start_US_sobol_anti_cv10;
	double price_US_sobol_anti_cv10 = USMC_sobol_anti_CV_Laguerre10.GetPrice();
	double var_US_sobol_anti_cv10 = USMC_sobol_anti_CV_Laguerre10.GetVariance();

	USMC_sobol_anti_CV_Laguerre10.OptimalNbSimul(tolerated_error);
	size_t Optimal_sobol_Nb_Simulation_anti_cv10 = USMC_sobol_anti_CV_Laguerre10.GetNbSimul();

	std::cout << " US Case - exec time for antithetic and control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynome: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_sobol_anti_cv10).count() << std::endl;
	std::cout << " US Case - price for antithetic and control variable MonteCarlo  with Sobol sequence and 6 order of Laguerre polynome : " << price_US_sobol_anti_cv10 << std::endl;
	std::cout << " US Case - variance for antithetic and control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynome : " << var_US_sobol_anti_cv10 << std::endl;
	std::cout << " US Case - optimal number of simulations forantithetic and control variable MonteCarlo with Sobol sequence and 6 order of Laguerre polynome : " << Optimal_sobol_Nb_Simulation_anti_cv10 << std::endl;

	delete sobol_gen_anti_cv10;
	delete ngen_sobol_anti_cv10;
	delete vectorG_sobol_anti_cv10;
	delete BSamerican_sobol_anti_cv10;
	return 0;

	}