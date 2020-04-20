#include "sde.hpp"

SinglePath::SinglePath(double start, double end, size_t nbSteps):
	StartTime(start),
	EndTime(end),
	Nbsteps(nbSteps)
{
	
};
void SinglePath::InsertValue(double val)
{
	Values.push_back(val);
};

double SinglePath::GetState(double time)
{
	size_t index = (size_t) (Nbsteps * (time - StartTime) / (EndTime - StartTime));
	return Values[index];
};

std::vector<double> SinglePath::GetAllValues()
{
	return Values;
}
///////////////////////////////////////////////////////////////////////////////////

RandomProcess::RandomProcess(RandomGenerator* Gen,int dim):
	Generator(Gen),
	dimension(dim)
{};

SinglePath* RandomProcess::GetPath(int pos)
{
	return Paths[pos];
}

////////////////////////////////////////////////////////////////////////////////////


Brownian1D::Brownian1D(RandomGenerator* Gen):
	RandomProcess(Gen,1)
{};

void Brownian1D::Simulate(double startTime,double endTime,size_t nbSteps)
{
	Paths.clear();
	SinglePath* Path = new SinglePath(startTime, endTime, nbSteps);
	Path->InsertValue(0.);
	double dt = (endTime - startTime) / nbSteps;
	double lastInserted = 0.;
	
	for(size_t i = 0; i < nbSteps; ++i)
	{
		double nextValue = lastInserted + sqrt(dt)* Generator->generate();
		Path->InsertValue(nextValue);
		lastInserted =nextValue;
	}
	Paths.push_back(Path);
	// delete Path;
};	
//////////////////////////////////////////////////////////////////////////////////////

BlackScholes1D::BlackScholes1D(RandomGenerator* Gen, double spot, double rate, double vol):
	Spot(spot),
	Rate(rate),
	Vol(vol),
	RandomProcess(Gen,1)
{};
///////////////////////////////////////////////////////////////////////////////////////
BSEuler1D::BSEuler1D(RandomGenerator* Gen, double spot, double rate, double vol):
	BlackScholes1D(Gen,spot,rate,vol)
{};

void BSEuler1D::Simulate(double startTime,double endTime,size_t nbSteps)
{
	Paths.clear();
	SinglePath* Path = new SinglePath(startTime, endTime, nbSteps);
	Path->InsertValue(Spot);
	double dt = (endTime - startTime) / nbSteps;
	double lastInserted = Spot;
	
	for(size_t i = 0; i < nbSteps; ++i)
	{
		//std::cout << i << std::endl;
		double nextValue = lastInserted + lastInserted * (Rate * dt + Vol * Generator->generate() * sqrt(dt));
		Path->InsertValue(nextValue);
		lastInserted =nextValue;
	}
	
	Paths.push_back(Path);
	// delete Path;
	
};
////////////////////////////////////////////////////////////////////////////////////////

BlackScholes2D::BlackScholes2D(RandomGenerator* Gen, double spot1,double spot2,double rate1, double rate2
		,double vol1,double vol2,double rho):
		Spot1(spot1),
		Rate1(rate1),
		Vol1(vol1),
		Spot2(spot2),
		Rate2(rate2),
		Vol2(vol2),
		Rho(rho),
		RandomProcess(Gen,2)
		
		
{};

BSEuler2D::BSEuler2D(RandomGenerator* Gen, double spot1,double spot2,double rate1, double rate2
		,double vol1,double vol2,double rho):
		BlackScholes2D(Gen, spot1, spot2, rate1, rate2
		,vol1,vol2,rho)
		
{};

void BSEuler2D::Simulate(double startTime,double endTime,size_t nbSteps)
{
	Paths.clear();
	
	SinglePath* Path1 = new SinglePath(startTime, endTime, nbSteps);
	SinglePath* Path2 = new SinglePath(startTime, endTime, nbSteps);
	
	Path1->InsertValue(Spot1);
	Path2->InsertValue(Spot2);
	
	double dt = (endTime - startTime) / nbSteps;
	double lastInserted1 = Spot1;
	double lastInserted2 = Spot2;

	for(size_t i = 0; i < nbSteps; ++i)
	{
		double E1 = Generator->generate();
		double Eind = Generator->generate();
		double E2 = Rho*E1 + sqrt(1-Rho*Rho)*Eind;
		
		double nextValue1 = lastInserted1 + lastInserted1 * (Rate1 * dt + Vol1 * E1 * sqrt(dt));
		double nextValue2 = lastInserted2 + lastInserted2 * (Rate2 * dt + Vol2 * E2 * sqrt(dt));
		Path1->InsertValue(nextValue1);
		Path2->InsertValue(nextValue2);
		lastInserted1 =nextValue1;
		lastInserted2 =nextValue2;
	}
	
	Paths.push_back(Path1);
	Paths.push_back(Path2);

	// delete Path;
		
};

// void export_csv(std::string f_name) const
// {
        // double sMin = s_mesh.get_Smin();
        // double dx = s_mesh.get_dx();
        
        // ofstream myFile;
		// myFile.open(f_name);
		
		// f << "Spot,Price,Delta,Gamma" << "\n";
		
		// for(int i = 0; i<solution.size(); ++i)
		// {
				// f <<  exp(sMin+i*dx) << "," << solution[i] << "," << delta[i] << "," << gamma[i] << "\n";
		// }

		// std::cout << "Results exported to " << f_name << std::endl;
		
        // f.close();
// }