#include "sde.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{		

/////////////////////////////////////////////////////////////////////////////////////
//INPUT FROM USER

	double startTime = 0.;
	double endTime = 1.;
	size_t nbsteps = 100;
	double dt = (endTime - startTime) / nbsteps;
//for 1D
	double spot = 100.;
	double vol= 0.2;
	double rate = 0.1;
//for 2D
	double spot1 = 100.;
	double spot2 = 100.;
	double vol1= 0.2;
	double vol2= 0.35;
	double rate1 = 0.01;	
	double rate2 = 0.03;
	double rho = -1.;
	
	UniformGenerator* ptr = new EcuyerCombined();	
	RandomGenerator* gtr = new NormalBoxMuller(ptr,0.,1.);
	
//for N Dim (3 assets to start)
std::vector<std::vector<double>> InitialSpot_vector = {{120},{100},{80}};

std::vector<std::vector<double>> Sigma_vector = {{0.25},{0.15},{0.3}};

std::vector<std::vector<double>> Mu_vector = {{0.05},{0.015},{0.005}};

std::vector<std::vector<double>> Correl_mat = {{1,-0.2,0.4},
											   {-0.2,1,0.2},
											   {0.4,0.2,1}};
											   
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
	matrix S(InitialSpot_vector);
	
	PayOffBasket* bsktcall = new PayOffBasketCall(W, S,100.);
	PayOffBasket* bsktcallCV = new PayOffControlVarBasketCall(W, S,100.);
	double payofftest = bsktcall->operator()(W,S);
	double payofftest2 = bsktcallCV->operator()(W,S);
	
	std::cout << "Payoff for the basket call is " << payofftest << std::endl;
	std::cout << "Payoff for the basket call Control Variate is " << payofftest2 << std::endl;
	delete bsktcall;
	delete bsktcallCV;
///////////////////////////////////////////////////////////////////////////////////////
//TEST MATRIX CLASS	

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

	BSEuler2D dynamics = BSEuler2D(gtr, spot1,spot2, rate1,rate2, vol1,vol2,rho);
	dynamics.Simulate(startTime, endTime, nbsteps);
	SinglePath* path1 = dynamics.GetPath(0);
	SinglePath* path2 = dynamics.GetPath(1);
	// std::cout << "debut print vector" << std::endl;
	// for(size_t i = 0; i< nbsteps;++i)
	// {
		// std::cout << "Step " << i << ", Spot= " << path1->GetState(i*dt) << ",";
	// }
	// std::cout << std::endl;
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Test GaussianVector generation

	matrix Sigma(Sigma_vector);
	matrix Nu(Mu_vector);
	matrix Correl(Correl_mat);

	matrix CovarMatrix = VarCovarMatrix(Sigma,Correl);
	
	Normal* ngnr = new NormalBoxMuller(ptr,0.,1.);
	
	// bool testinvertible = isInvertible(CovarMatrix,CovarMatrix.nb_rows());
	// std::cout<< "test invertivle returns: " << testinvertible << std::endl;
	GaussianVector* gvec;
	
	// if(testinvertible == 1)
	gvec = new GaussianVectorCholesky(ngnr,Nu,Sigma,Correl,CovarMatrix);
	// else
		// gvec = new GaussianVectorDiag(ngnr,Nu,Sigma,Correl,CovarMatrix);

	matrix output = gvec->CorrelatedGaussianVector();
	std::cout << "Output" <<std::endl;
	output.Print();
	matrix output2 = gvec->CorrelatedGaussianVector();
	std::cout << "Output2" <<std::endl;
	output2.Print();
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
	delete path1;
	delete path2;
	
	delete ngnr;
	delete gvec;

	return 0;
}