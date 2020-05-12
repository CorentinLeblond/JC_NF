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
	double dt = (endTime - startTime) / nbsteps;

	size_t N = 1000;
/*
//for 1D
	double spot = 100.;
	double vol= 0.2;
	double rate = 0.05;
//for 2D
	double spot1 = 100.;
	double spot2 = 100.;
	double vol1= 0.2;
	double vol2= 0.2;
	double rate1 = 0.00;	
	double rate2 = 0.03;
	double rho = -1.;
	
	UniformGenerator* ptr = new EcuyerCombined();	
	RandomGenerator* gtr = new NormalBoxMuller(ptr,0.,1.);
*/	
//for N Dim (3 assets to start)
double rate = 0.05;
std::vector<std::vector<double>> Spot_vector = {{100},{120},{80}};

std::vector<std::vector<double>> Sigma_vector = {{0.25},{0.20},{0.15}};

std::vector<std::vector<double>> Mu_vector = {{rate},{rate},{rate}};

std::vector<std::vector<double>> Correl_mat = {{1,-0.2,0.4},
											   {-0.2,1,0.4},
											   {0.4,0.4,1}};
											   
std::vector<std::vector<double>> Weights_mat ={{0.3,0.5,0.2}};
//////////////////////////////////////////////////////////////////////////////////////
//Test Van der Corput
	//UniformGenerator* qsr = new VanDerCorput(2,1);
	//
	//double test = 0.;
	//
	//for(size_t i = 0;i<10;++i)
	//{
	//	test = qsr->generate();
	//	std::cout << "Generated quasi rdm nbr is: " << test << std::endl;
	//}
	//delete qsr;
///////////////////////////////////////////////////////////////////////////////////////
//Test payoff
	matrix W(Weights_mat);
	matrix S(Spot_vector);
	
	PayOffBasket* bsktcall = new PayOffBasketCall(W, S,100.);
	PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(W, S,100.);
	double payofftest = bsktcall->operator()(S);
	double payofftest2 = bsktcallCV->operator()(S);
	
	std::cout << "Payoff for the basket call is " << payofftest << std::endl;
	std::cout << "Payoff for the basket call Control Variate is " << payofftest2 << std::endl;

////////////////////////////////////////////////////////////////////////////////////////
	
//First thing we want to do is to test that our EcuyerCombined method really generates a uniform distribution (0,1)
	
	// double mean_unif = ptr -> mean(100000);
	// double var_unif = ptr -> var(100000);
	
	// std::cout << "Mean Uniform : " << mean_unif << std::endl;
	// std::cout << "Var Uniform : " << var_unif << std::endl;
	
	// std::ofstream myFile1;
	// myFile1.open("uniform.csv");
	
	// for(int i = 0; i< 10000;++i)
	// {
		// myFile1 << ptr->generate() << std::endl;
	// }
	
//Second thing we want to test that our Box muller is fine
	
	// double mean_normal = gtr -> mean(100000);
	// double var_normal = gtr -> var(100000);
	// std::cout << "Mean Normal : " << mean_normal << std::endl;
	// std::cout << "Var Normal : " << var_normal << std::endl;
	// std::ofstream myFile2;
	// myFile2.open("Normal.csv");
	
	// for(int i = 0; i< 10000;++i)
	// {
		// myFile2 << gtr->generate() << std::endl;
	// }
	
	//RandomProcess* rdm = new Brownian1D(gtr);
///////////////////////////////////////////////////////////////////////////////////////////////	
//Test BS1D

	// BSEuler1D dynamics = BSEuler1D(gtr, spot, rate, vol);
	// dynamics.Simulate(startTime, endTime, nbsteps);
	// SinglePath* path = dynamics.GetPath(0);
///////////////////////////////////////////////////////////////////////////////////////////////	
//Test BS2D

	//BSEuler2D dynamics = BSEuler2D(gtr, spot1,spot2, rate1,rate2, vol1,vol2,rho);
	//dynamics.Simulate(startTime, endTime, nbsteps);
	//SinglePath* path1 = dynamics.GetPath(0);
	//SinglePath* path2 = dynamics.GetPath(1);

///////////////////////////////////////////////////////////////////////////////////////////////
//Test BSND

	matrix spot_m(Spot_vector);
	matrix Sigma(Sigma_vector);
	matrix Nu(Mu_vector);
	matrix Nu2(Mu_vector);
	matrix Correl(Correl_mat);



	Nu*= dt;

	// std::cout << " Nu after dt operator " << std::endl;
	// Nu.Print();

	matrix CovarMatrix = VarCovarMatrix(Sigma, Correl);

	double K = 100;

	ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(W, spot_m, CovarMatrix, K,
		rate, endTime);

	//std::vector<std::vector<double>> Spot_vector_maturity = { {130},{110},{91} };
	//matrix endspot(Spot_vector_maturity);
	double df = exp(-rate * endTime);
	double price_cf = CFbaskt->operator()(spot_m, df);
	std::cout << "closed formula for the basket option " <<  price_cf << std::endl;
	
	UniformGenerator* ugen = new EcuyerCombined();
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);
	GaussianVectorCholesky* corrG = new GaussianVectorCholesky(ngen, Sigma, Correl, CovarMatrix);
	RandomProcess* path = new BSEulerND(corrG,spot_m,rate);

	
	
	clock::time_point start = clock::now(); //We start the chrono at that point of the code
	EuropeanBasket MC(N, bsktcall, path);
	MC.Simulate(startTime,endTime,nbsteps);
	clock::time_point end = clock::now(); //We take again the time once the entire simulation is done
	clock::duration execution_time = end - start; //We compute the differentce and print it next line
	//Format is in seconds (cast <std::ratio<1>) stands for seconds
	std::cout << "exec time for vanilla Basket call: " << std::chrono::duration <double,std::ratio<1>> (execution_time).count() << std::endl;
	//Les variables de temps sont déclarées, on peut les réutiliser plus loin dans le code comme n'importe quelle autre variable

	double priceMC = MC.GetPrice();
	
	//dynamics.Simulate(startTime, endTime, nbsteps);
	//matrix chemin = dynamics.GetAllPaths();
	//chemin.Print();

	UniformGenerator* ugen2 = new EcuyerCombined();
	RandomGenerator* ngen2 = new NormalBoxMuller(ugen2, 0., 1.);
	// RandomProcess* chemin = new BSEuler1D(ngen2, spot, rate, vol);
/*
	

	PayOffCall* call = new PayOffCall(K);

	EuropeanVanilla_MonteCarlo VMC(N, call, chemin);
	VMC.Simulate(startTime, endTime, nbsteps);
	double priceVMC = VMC.GetPrice(rate, endTime);
*/	

	UniformGenerator* ugen3 = new EcuyerCombined();
	Normal* ngen3 = new NormalBoxMuller(ugen3, 0., 1.);
	GaussianVectorCholesky* corrGauss = new GaussianVectorCholesky(ngen3, Sigma, Correl, CovarMatrix);
	RandomProcess* path_cv = new BSEulerND(corrGauss,spot_m,rate);

	clock::time_point start_ = clock::now();
	EuropeanBasket_controlvariable CVMC(N, bsktcall, bsktcallCV, path_cv,price_cf);
	CVMC.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_ = clock::now(); //We take again the time once the entire simulation is done
	clock::duration execution_time2 = end_ - start_; //We compute the differentce and print it next line
	//Format is in seconds (cast <std::ratio<1>) stands for seconds
	std::cout << "exec time for controle variate Basket call: " << std::chrono::duration <double, std::ratio<1>>(execution_time2).count() << std::endl;
	double priceCVMC = CVMC.GetPrice();

	// double varVMC = VMC.GetVariance();
	double varMC = MC.GetVariance();
	double varCVMC = CVMC.GetVariance();

	// std::cout << "price Vanilla 1D MC " << priceVMC << std::endl;
	std::cout << "price Vanilla CV MC " << priceCVMC << std::endl;
	std::cout << "price Vanilla MC for Basket " << priceMC << std::endl;

	// std::cout << "variance Vanilla 1D MC " << varVMC << std::endl;
	std::cout << "variance Vanilla CV MC " << varCVMC << std::endl;
	std::cout << "variance Vanilla MC for Basket " << varMC << std::endl;

	UniformGenerator* vdc_gen = new VanDerCorput(2, 1);
	Normal* ngen4 = new NormalBoxMuller(vdc_gen, 0., 1.);
	GaussianVectorCholesky* GaussVDC = new GaussianVectorCholesky(ngen4, Sigma, Correl, CovarMatrix);
	RandomProcess* BS_vdc_vc = new BSEulerND(GaussVDC, spot_m, rate);

	EuropeanBasket_controlvariable MC_quasi_vc(N, bsktcall, bsktcallCV, BS_vdc_vc,price_cf);
	MC_quasi_vc.Simulate(startTime, endTime, nbsteps);
	double price_quasi_vc = MC_quasi_vc.GetPrice();
	double var_quasi = MC_quasi_vc.GetVariance();

	std::cout << "Price of quasi and VC " << price_quasi_vc << std::endl;
	std::cout << "Variance quasi and VC " << var_quasi << std::endl;


	UniformGenerator* ugen5 = new EcuyerCombined();
	Normal* ngen5 = new NormalBoxMuller(ugen5, 0., 1.);
	GaussianVectorCholesky* Gauss_anti = new GaussianVectorCholesky(ngen5, Sigma, Correl, CovarMatrix);
	RandomProcess* BS_anti = new BSEulerNDAntithetic(Gauss_anti, spot_m, rate);

	EuropeanBasket_Antithetic MC_anti(N, bsktcall, BS_anti);
	MC_anti.Simulate(startTime, endTime, nbsteps);
	double price_anti = MC_anti.GetPrice();
	double var_anti = MC_anti.GetVariance();

	std::cout << "Price of anti " << price_anti<< std::endl;
	std::cout << "Variance anti " << var_anti << std::endl;

//////// TEST US ////////////////////////////////////////////////////////////

	std::vector<basis_functions*> basefunc;
	
	basefunc.push_back(new Poly_Laguerre(0));
	basefunc.push_back(new Poly_Laguerre(1));
	basefunc.push_back(new Poly_Laguerre(2));
	//basefunc.push_back(new Laguerre_test(0));
	//basefunc.push_back(new Laguerre_test(1));
	//basefunc.push_back(new Laguerre_test(2));


	double df2 = exp(-rate * dt);

	UniformGenerator* ugen6 = new EcuyerCombined();
	Normal* ngen6 = new NormalBoxMuller(ugen6, 0., 1.);

	GaussianVectorCholesky* vectorG = new GaussianVectorCholesky(ngen6, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican = new BSEulerND(vectorG, spot_m, rate);

	clock::time_point start_US = clock::now();
	AmericanMonteCarlo_basket USMC(N, bsktcall,BSamerican,basefunc);
	USMC.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US = clock::now(); //We take again the time once the entire simulation is done
	clock::duration execution_timeUS = end_US - start_US; //We compute the differentce and print it next line
	std::cout << "exec time for US MC Basket call: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS).count() << std::endl;
	double price_US = USMC.GetPrice();
	double var_US = USMC.GetVariance();
	//matrix early = USMC.GetEarlyExec();

	//std::cout << "Early exercice " << price_US << std::endl;
	//early.Print();
	std::cout << "Price of US " << price_US << std::endl;
	std::cout << "Variance US " << var_US << std::endl;


	UniformGenerator* ugen7 = new EcuyerCombined();
	Normal* ngen7 = new NormalBoxMuller(ugen7, 0., 1.);

	GaussianVectorCholesky* vectorG7 = new GaussianVectorCholesky(ngen7, Sigma, Correl, CovarMatrix);
	RandomProcess* BSamerican_CV = new BSEulerND(vectorG7, spot_m, rate);

	clock::time_point start_US_CV = clock::now();
	AmericanMonteCarlo_basket_controlevariable USMC_CV(N, bsktcall, bsktcallCV, BSamerican_CV, basefunc,price_cf);
	USMC_CV.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_CV = clock::now(); //We take again the time once the entire simulation is done
	clock::duration execution_timeUS_CV = end_US_CV - start_US_CV; //We compute the differentce and print it next line
	std::cout << "exec time for US MC Basket call with control variate: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_CV).count() << std::endl;
	double price_US_CV = USMC_CV.GetPrice();
	double var_US_CV = USMC_CV.GetVariance();



	UniformGenerator* ugen8 = new EcuyerCombined();
	Normal* ngen8 = new NormalBoxMuller(ugen8, 0., 1.);

	GaussianVectorCholesky* vectorG8 = new GaussianVectorCholesky(ngen8, Sigma, Correl, CovarMatrix);
	RandomProcess* BS_anti_US = new BSEulerNDAntithetic(vectorG8, spot_m, rate);

	clock::time_point start_US_anti = clock::now();
	AmericanMonteCarlo_basket_Antithetic USMC_anti(N, bsktcall, BS_anti_US, basefunc);
	USMC_anti.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_anti = clock::now(); //We take again the time once the entire simulation is done
	clock::duration execution_timeUS_anti = end_US_anti - start_US_anti; //We compute the differentce and print it next line
	std::cout << "exec time for US MC Basket call with antitethic: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti).count() << std::endl;
	double price_US_anti = USMC_anti.GetPrice();
	double var_US_anti = USMC_anti.GetVariance();
	//matrix early = USMC.GetEarlyExec();


	UniformGenerator* ugen9 = new EcuyerCombined();
	Normal* ngen9 = new NormalBoxMuller(ugen9, 0., 1.);

	GaussianVectorCholesky* vectorG9 = new GaussianVectorCholesky(ngen9, Sigma, Correl, CovarMatrix);
	RandomProcess* BS_anti_CV_US = new BSEulerNDAntithetic(vectorG9, spot_m, rate);

	clock::time_point start_US_CV_anti = clock::now();
	AmericanMonteCarlo_basket_Antithetic_CV USMC_anti_CV(N, bsktcall, bsktcallCV, BS_anti_CV_US, basefunc, price_cf);
	USMC_anti_CV.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_US_anti_CV = clock::now(); //We take again the time once the entire simulation is done
	clock::duration execution_timeUS_anti_CV = end_US_anti_CV - start_US_CV_anti; //We compute the differentce and print it next line
	std::cout << "exec time for US MC Basket call with CV and antitethic: " << std::chrono::duration <double, std::ratio<1>>(execution_timeUS_anti_CV).count() << std::endl;
	double price_US_anti_CV = USMC_anti_CV.GetPrice();
	double var_US_anti_CV = USMC_anti_CV.GetVariance();
	//matrix early = USMC.GetEarlyExec();



	std::cout << "Price of US CV " << price_US_CV << std::endl;
	std::cout << "Variance US CV " << var_US_CV << std::endl;

	std::cout << "Price of US anti " << price_US_anti << std::endl;
	std::cout << "Variance US anti " << var_US_anti << std::endl;

	std::cout << "Price of US CV anti " << price_US_anti_CV << std::endl;
	std::cout << "Variance US CV anti " << var_US_anti_CV << std::endl;



	//////////////// BERMUDEAN TEST ////////////////////////////////////////

	//std::vector<basis_functions*> basefunc;
	//
	//basefunc.push_back(new Poly_Laguerre(1));
	//basefunc.push_back(new Poly_Laguerre(2));

	matrix Schedule_exec(10, 1);

	double test = 0.1;

	for (size_t i = 0; i < 10; i++) 
	{
		
		Schedule_exec(i, 0) = test;
		test = test + 0.05;
	}

	//std::cout << "Exec possibles" << std::endl;

	//Schedule_exec.Print();
	size_t wk = 0;
	CalendarManagement* wkdayfunc = new rounded_workingdays(wk);

	UniformGenerator* ugen10 = new EcuyerCombined();
	Normal* ngen10 = new NormalBoxMuller(ugen10, 0., 1.);

	GaussianVectorCholesky* vectorG10 = new GaussianVectorCholesky(ngen10, Sigma, Correl, CovarMatrix);
	RandomProcess* BSberm = new BSEulerND(vectorG10, spot_m, rate);

	clock::time_point start_berm = clock::now();
	Bermudean_BasketOption berm(N, bsktcall, BSberm, basefunc, wkdayfunc, Schedule_exec);
	berm.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_berm = clock::now(); //We take again the time once the entire simulation is done
	clock::duration execution_timeberm = end_berm - start_berm; //We compute the differentce and print it next line
	std::cout << "exec time for berm MC Basket call: " << std::chrono::duration <double, std::ratio<1>>(execution_timeberm).count() << std::endl;
	double price_berm = berm.GetPrice();
	double var_berm = berm.GetVariance();
	//matrix early = USMC.GetEarlyExec();

	//std::cout << "Early exercice " << price_US << std::endl;
	//early.Print();
	std::cout << "Price of berm " << price_berm << std::endl;
	std::cout << "Variance berm " << var_berm << std::endl;


	UniformGenerator* ugen11 = new EcuyerCombined();
	Normal* ngen11 = new NormalBoxMuller(ugen11, 0., 1.);

	GaussianVectorCholesky* vectorG11 = new GaussianVectorCholesky(ngen11, Sigma, Correl, CovarMatrix);
	RandomProcess* BSberm_CV = new BSEulerND(vectorG11, spot_m, rate);

	clock::time_point start_berm_CV = clock::now();
	Bermudean_BasketOption_CV berm_CV(N, bsktcall, bsktcallCV, BSberm_CV, basefunc, wkdayfunc, Schedule_exec, price_cf);
	berm_CV.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_berm_CV = clock::now(); //We take again the time once the entire simulation is done
	clock::duration execution_timebermCV = end_berm_CV - start_berm_CV; //We compute the differentce and print it next line
	std::cout << "exec time for berm CV MC Basket call: " << std::chrono::duration <double, std::ratio<1>>(execution_timebermCV).count() << std::endl;
	double price_berm_CV = berm_CV.GetPrice();
	double var_berm_CV = berm_CV.GetVariance();

	std::cout << "Price of berm CV " << price_berm_CV << std::endl;
	std::cout << "Variance berm CV " << var_berm_CV << std::endl;



	UniformGenerator* ugen12 = new EcuyerCombined();
	Normal* ngen12 = new NormalBoxMuller(ugen12, 0., 1.);

	GaussianVectorCholesky* vectorG12 = new GaussianVectorCholesky(ngen12, Sigma, Correl, CovarMatrix);
	RandomProcess* BS_berm_anti = new BSEulerNDAntithetic(vectorG12, spot_m, rate);

	clock::time_point start_berm_anti = clock::now();
	Bermudean_BasketOption_antithetic berm_anti(N, bsktcall, BS_berm_anti, basefunc, wkdayfunc, Schedule_exec);
	berm_anti.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_berm_anti = clock::now(); //We take again the time once the entire simulation is done
	clock::duration execution_timebermanti = end_berm_anti - start_berm_anti; //We compute the differentce and print it next line
	std::cout << "exec time for berm anti MC Basket call: " << std::chrono::duration <double, std::ratio<1>>(execution_timebermanti).count() << std::endl;
	double price_berm_anti = berm_anti.GetPrice();
	double var_berm_anti = berm_anti.GetVariance();

	std::cout << "Price of berm anti " << price_berm_anti << std::endl;
	std::cout << "Variance berm anti " << var_berm_anti << std::endl;


	UniformGenerator* ugen13 = new EcuyerCombined();
	Normal* ngen13 = new NormalBoxMuller(ugen13, 0., 1.);

	GaussianVectorCholesky* vectorG13 = new GaussianVectorCholesky(ngen13, Sigma, Correl, CovarMatrix);
	RandomProcess* BS_berm_anti_CV = new BSEulerNDAntithetic(vectorG13, spot_m, rate);

	clock::time_point start_berm_anti_CV = clock::now();
	Bermudean_BasketOption_antithetic_CV berm_anti_CV(N, bsktcall, bsktcallCV, BS_berm_anti_CV, basefunc, wkdayfunc, Schedule_exec, price_cf);
	berm_anti_CV.Simulate(startTime, endTime, nbsteps);
	clock::time_point end_berm_anti_CV = clock::now(); //We take again the time once the entire simulation is done
	clock::duration execution_timebermantiCV = end_berm_anti_CV - start_berm_anti_CV; //We compute the differentce and print it next line
	std::cout << "exec time for berm anti MC Basket call: " << std::chrono::duration <double, std::ratio<1>>(execution_timebermantiCV).count() << std::endl;
	double price_berm_anti_CV = berm_anti_CV.GetPrice();
	double var_berm_anti_CV = berm_anti_CV.GetVariance();

	std::cout << "Price of berm anti and CV " << price_berm_anti_CV << std::endl;
	std::cout << "Variance berm anti and CV " << var_berm_anti_CV << std::endl;


	
	//std::vector<std::vector<double>> vecteur_T = { {3.,-339.878441265584,18958.7374809538},
	//{-339.878441265584,38600.2318444388,-2158526.16480462},
	//{18958.7374809538,-2158526.16480462,121009687.632087} };

	//std::vector<std::vector<double>> vecteur = { {1.,2.},{4.,5.},{7.,15.} };

	//matrix MM(vecteur_T);

	//matrix I = Inverse(MM, MM.nb_rows());
	////matrix C = Inverse_Cholesky(vecteur);
	//I.Print();
	//std::cout << I(0, 0) << std::endl;

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

	// delete ptr;
	// delete gtr;
	//delete corrG;
	//delete corrGauss;
	//delete path1;
	//delete path2;

	////delete vdc_gen;
	delete ngen4;
	delete ugen3;
	delete ugen5;
	delete ngen3;
	delete ngen5;
	delete ugen6;
	delete ngen6;
	delete ugen7;
	delete ngen7;
	delete ugen8;
	delete ngen8;
	delete ugen9;
	delete ngen9;
	delete ugen10;
	delete ngen10;
	delete ugen11;
	delete ngen11;
	delete ugen12;
	delete ngen12;
	delete ugen13;
	delete ngen13;
	for (size_t i = 0; i < basefunc.size(); i++) {
		delete basefunc[i];
	};

	delete path;
	delete ugen;
	delete ngen;
	delete ugen2;
	delete ngen2;

	delete bsktcall;
	delete bsktcallCV;
	
	//delete ngnr;
	//delete gvec;

	return 0;
}