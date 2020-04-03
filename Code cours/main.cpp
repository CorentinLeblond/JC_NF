#include "sde.hpp"
#include <iostream>
#include <fstream>

int main(int argc, char* argv[])
{		

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
		// std::cout << "Step " << i << ", Spot= " << path->GetState(i*dt) << ",";
	// }
	// std::cout << std::endl;
	
	std::ofstream myFile;
	myFile.open("BlackScholes2D_3.csv");
	
	for(int i = 0; i< nbsteps;++i)
	{
		myFile << i*dt << "," << path1->GetState(i*dt) << "," << path2->GetState(i*dt)<<std::endl;
	}
	
	delete ptr;
	delete gtr;
	delete path1;
	delete path2;
	
	return 0;
}