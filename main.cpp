#include "sde.hpp"
#include "MonteCarlo.h"
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{		

/////////////////////////////////////////////////////////////////////////////////////
//INPUT FROM USER

	double startTime = 0.;
	double endTime = 1.;
	size_t nbsteps = 252;
	double dt = (endTime - startTime) / nbsteps;
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
	
//for N Dim (3 assets to start)
std::vector<std::vector<double>> Spot_vector = {{100},{120},{80}};

std::vector<std::vector<double>> Sigma_vector = {{0.25},{0.20},{0.15}};

std::vector<std::vector<double>> Mu_vector = {{rate},{rate},{rate}};

std::vector<std::vector<double>> Correl_mat = {{1,-0.2,0.4},
											   {-0.2,1,0.4},
											   {0.4,0.4,1}};
											   
std::vector<std::vector<double>> Weights_mat ={{0.3,0.5,0.2}};
//////////////////////////////////////////////////////////////////////////////////////
//Test Van der Corput
	UniformGenerator* qsr = new VanDerCorput(2,1);
	
	double test = 0.;
	
	for(size_t i = 0;i<10;++i)
	{
		test = qsr->generate();
		std::cout << "Generated quasi rdm nbr is: " << test << std::endl;
	}
	delete qsr;
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

///////////////////////////////////////////////////////////////////////////////////////
//TEST MATRIX CLASS	
		//matrix y({{120},{100},{80}});
		//double test_mean = y.mean();
		//double test_variance = y.variance();
		
		//std::cout << "Mean is " << test_mean << std::endl;
		//std::cout << "Variance is " << test_variance << std::endl;
		// matrix y(3,3);
		// y(0,0) = 4.;
		// y(1,0)  = 12.;
		// y(2,0) = -16.;
		// y(0,1)  = 12.;		
		// y(1,1) = 37;
		// y(1,2)  = -43;
		// y(2,0)  = -16.;		
		// y(2,1) = -43;
		// y(2,2)  = 98;		
		
		// std::cout << "y matrix" << std::endl;
		// y.Print();
		// matrix o = y.Cholesky();
		// std::cout << "o matrix" << std::endl;
		// o.Print();
		
		// matrix x(2,1);
		// x(0,0) = 1.;
		// x(1,0)  = 2.;
		// x.Print();
		// x.Resize(2,2);
		
		// x.Diagonalization();
		// x.Print();
		
		// matrix m(2, 2);
		// matrix r(2, 2);
        // // std::cout << m.nb_rows() << std::endl;
        // // std::cout << m.nb_cols() << std::endl;
		
		// m(0,0) = 1.;
		// m(0,1)  = 2.;
		// m(1,0)  = 3.;
		// m(1,1)  = 4.;
		// std::cout << "matrix m" << std::endl;
		// m.Print();
		// std::cout << "matrix r" << std::endl;
		// r(0,0) = 1.;
		// r(0,1)  = 2.;
		// r(1,0)  = 3.;
		// r(1,1)  = 4.;
		// r.Print();
		// std::cout << "matrix z" << std::endl;
		// matrix z = m * r;
		// z.Print();
		
		// std::cout << "matrix e" << std::endl;
		// // matrix e = m * r;
		// // e.Print();
		// // m+=r;
		// // m.Print();
		
		// // matrix g = m + r;
		// // g.Print();
		// std::cout << "matrix m" << std::endl;
		// m*=r;
		// m.Print();
		// // m.Print();
		
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////		
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
	
//Test BS2D

	//BSEuler2D dynamics = BSEuler2D(gtr, spot1,spot2, rate1,rate2, vol1,vol2,rho);
	//dynamics.Simulate(startTime, endTime, nbsteps);
	//SinglePath* path1 = dynamics.GetPath(0);
	//SinglePath* path2 = dynamics.GetPath(1);


//Test BSND

	UniformGenerator* ugen = new EcuyerCombined();
	Normal* ngen = new NormalBoxMuller(ugen, 0., 1.);

	matrix spot_m(Spot_vector);
	matrix Sigma(Sigma_vector);
	matrix Nu(Mu_vector);
	matrix Nu2(Mu_vector);
	matrix Correl(Correl_mat);

	Nu*= dt;

	std::cout << " Nu after dt operator " << std::endl;
	Nu.Print();

	matrix CovarMatrix = VarCovarMatrix(Sigma, Correl);

	GaussianVectorCholesky* corrG = new GaussianVectorCholesky(ngen, Nu, Sigma, Correl, CovarMatrix);

	RandomProcess* path = new BSEulerND(corrG,spot_m);

	size_t N = 10000;

	EuropeanBasket_MonteCarlo MC(N, bsktcall, path);
	MC.Simulate(startTime,endTime,nbsteps);
	double priceMC = MC.GetPrice(rate, endTime);
	//dynamics.Simulate(startTime, endTime, nbsteps);
	//matrix chemin = dynamics.GetAllPaths();
	//chemin.Print();

	UniformGenerator* ugen2 = new EcuyerCombined();
	RandomGenerator* ngen2 = new NormalBoxMuller(ugen2, 0., 1.);

	RandomProcess* chemin = new BSEuler1D(ngen2, spot, rate, vol);

	double K = 100;

	PayOffCall* call = new PayOffCall(K);

	EuropeanVanilla_MonteCarlo VMC(N, call, chemin);
	VMC.Simulate(startTime, endTime, nbsteps);
	double priceVMC = VMC.GetPrice(rate, endTime);
	

	UniformGenerator* ugen3 = new EcuyerCombined();
	Normal* ngen3 = new NormalBoxMuller(ugen3, 0., 1.);

	GaussianVectorCholesky* corrGauss = new GaussianVectorCholesky(ngen3, Nu, Sigma, Correl, CovarMatrix);

	RandomProcess* path_cv = new BSEulerND(corrGauss,spot_m);

	EuropeanBasket_MonteCarlo_controlvariable CVMC(N, bsktcall, bsktcallCV, path_cv);
	CVMC.Simulate(startTime, endTime, nbsteps);
	double priceCVMC = CVMC.GetPrice(rate, endTime);

	double varVMC = VMC.GetVariance();
	double varMC = MC.GetVariance();
	double varCVMC = CVMC.GetVariance();

	std::cout << "price Vanilla 1D MC " << priceVMC << std::endl;
	std::cout << "price Vanilla CV MC " << priceCVMC << std::endl;
	std::cout << "price Vanilla MC for Basket " << priceMC << std::endl;

	std::cout << "variance Vanilla 1D MC " << varVMC << std::endl;
	std::cout << "variance Vanilla CV MC " << varCVMC << std::endl;
	std::cout << "variance Vanilla MC for Basket " << varMC << std::endl;
	
	ClosedFormulaBasketCall* CFbaskt = new ClosedFormulaBasketCall(W, spot_m, CovarMatrix, K,
		rate, endTime);

	//std::vector<std::vector<double>> Spot_vector_maturity = { {130},{110},{91} };
	//matrix endspot(Spot_vector_maturity);
	double df = exp(-rate * endTime);
	std::cout << "closed formula for the basket option" << CFbaskt->operator()(spot_m,df) << std::endl;
	
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
////////////////////////////////////////////////////////////////////////////////////////	
//TEST CSV

	// std::ofstream myFile;
	// myFile.open("BlackScholes2D_3.csv");
	
	// for(int i = 0; i< nbsteps;++i)
	// {
		// myFile << i*dt << "," << path1->GetState(i*dt) << "," << path2->GetState(i*dt)<<std::endl;
	// }
	
	delete ptr;
	delete gtr;
	delete corrG;
	delete corrGauss;
	//delete path1;
	//delete path2;

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