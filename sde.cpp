#include "sde.hpp"

///////////////////////////////SINGLE PATH CLASS //////////////////////////////////////////
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

///////////////////////////////RANDOM PROCESS CLASS //////////////////////////////////////////

RandomProcess::RandomProcess(RandomGenerator* Gen,int dim):
	Generator(Gen),
	dimension(dim)
{};

SinglePath* RandomProcess::GetPath(int pos)
{
	return Paths[pos];
}

matrix RandomProcess::GetAllPaths() {


	std::vector<double> element;
	std::vector<std::vector<double>> res;
	for (size_t single = 0; single < Paths.size(); single++) {

		element = Paths[single]->GetAllValues();
		res.push_back(element);

		}

	matrix chemins(res);
	return chemins;

};


matrix RandomProcess::GetAllPathsAnti()
{
	std::vector<double> element;
	std::vector<std::vector<double>> res;
	for (size_t single = 0; single < PathsAntithetic.size(); single++)
	{

		element = PathsAntithetic[single]->GetAllValues();
		res.push_back(element);

	}

	matrix chemins(res);
	return chemins;

};

double RandomProcess::Get_Dt()
{
	return dt_sde;
};

double RandomProcess::Get_rate()
{
	return rate_sde;
};


Brownian1D::Brownian1D(RandomGenerator* Gen):
	RandomProcess(Gen,1)
{};

void Brownian1D::Simulate(double startTime,double endTime,size_t nbSteps)
{
	Paths.clear();
	SinglePath* Path = new SinglePath(startTime, endTime, nbSteps);
	Path->InsertValue(0.);
	double dt = (endTime - startTime) / nbSteps;
	dt_sde = dt;
	double lastInserted = 0.;
	
	for(size_t i = 0; i < nbSteps; ++i)
	{
		double nextValue = lastInserted + sqrt(dt)* Generator->generate();
		Path->InsertValue(nextValue);
		lastInserted =nextValue;
	}
	Paths.push_back(Path);
};	


BlackScholes1D::BlackScholes1D(RandomGenerator* Gen, double spot, double rate, double vol):
	Spot(spot),
	Rate(rate),
	Vol(vol),
	RandomProcess(Gen,1)
{
	rate_sde = rate;

};

BSEuler1D::BSEuler1D(RandomGenerator* Gen, double spot, double rate, double vol):
	BlackScholes1D(Gen,spot,rate,vol)
{};

void BSEuler1D::Simulate(double startTime,double endTime,size_t nbSteps)
{
	Paths.clear();
	SinglePath* Path = new SinglePath(startTime, endTime, nbSteps);
	Path->InsertValue(Spot);
	double dt = (endTime - startTime) / nbSteps;
	dt_sde = dt;
	double lastInserted = Spot;
	
	for(size_t i = 0; i < nbSteps; ++i)
	{
		double nextValue = lastInserted + lastInserted * (Rate * dt + Vol * Generator->generate() * sqrt(dt));
		Path->InsertValue(nextValue);
		lastInserted =nextValue;
	}
	
	Paths.push_back(Path);
	
};

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
	dt_sde = dt;
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
		
};

BlackScholesND::BlackScholesND(GaussianVector* CorrelGaussian, matrix spot_vec,double inputrate)
	:m_gaussian(CorrelGaussian),
	V_spot(spot_vec),
	rate(inputrate)
{
	rate_sde = inputrate;
};

BSEulerND::BSEulerND(GaussianVector* CorrelGaussian, matrix spot_vec,double inputrate):
	BlackScholesND(CorrelGaussian,spot_vec,inputrate)
{ 
	rate_sde = inputrate;
};


void BSEulerND::Simulate(double startTime, double endTime, size_t nbSteps)
{
	Paths.clear();
	size_t assets = V_spot.nb_rows();
	double last_spot = 0.;
	double next_spot = 0.;
	double dt = (endTime - startTime) / nbSteps;
	dt_sde = dt;
	Brownian.Resize(assets, nbSteps);
	matrix mean_vector(assets, 1);

	for (size_t t = 0; t < nbSteps; ++t)
	{
		matrix X = m_gaussian->CorrelatedGaussianVector();

		for (size_t i = 0; i < assets; ++i)
		{
			//create the matrix of all brownian motions 
			Brownian(i, t) = X(i, 0);
		}
	}

	Paths.resize(assets);

	for (size_t i = 0; i < assets; ++i)
	{
		Paths[i] = new SinglePath(startTime, endTime, nbSteps);
		
		last_spot = V_spot(i, 0);

		next_spot = 0.;
		
		Paths[i]->InsertValue(last_spot);
		
		for (size_t t = 0; t < nbSteps; ++t)
		{
			next_spot = last_spot * (1+ rate*dt + sqrt(dt)*Brownian(i, t));
			Paths[i]->InsertValue(next_spot);
			last_spot = next_spot;

		}

	}

};

BSEulerNDAntithetic::BSEulerNDAntithetic(GaussianVector* CorrelGaussian, matrix spot_vec,double inputrate):
	BlackScholesND(CorrelGaussian,spot_vec,inputrate)
{ 
	rate_sde = inputrate;
};


void BSEulerNDAntithetic::Simulate(double startTime, double endTime, size_t nbSteps)
{
	Paths.clear();
	PathsAntithetic.clear();
	
	size_t assets = V_spot.nb_rows();
	double last_spot = 0.;
	double last_spotAnti = 0.;
	double next_spot = 0.;
	double next_spotAnti = 0.;
	double dt = (endTime - startTime) / nbSteps;
	dt_sde = dt;
	
	Brownian.Resize(assets, nbSteps);
	BrownianAntithetic.Resize(assets, nbSteps);

	for (size_t t = 0; t < nbSteps; ++t)
	{
		matrix X = m_gaussian->CorrelatedGaussianVector();
		matrix XAntithetic = X*(-1);
		
		for (size_t i = 0; i < assets; ++i)
		{
			//create the matrix of all brownian motion 
			Brownian(i, t) = X(i, 0);
			BrownianAntithetic(i, t) = XAntithetic(i, 0);
		}

	}

	Paths.resize(assets);
	PathsAntithetic.resize(assets);
	
	for (size_t i = 0; i < assets; ++i)
	{

		Paths[i] = new SinglePath(startTime, endTime, nbSteps);
		PathsAntithetic[i] = new SinglePath(startTime, endTime, nbSteps);
	
		last_spot = V_spot(i, 0);
		last_spotAnti = V_spot(i, 0);

		next_spot = 0.;
		next_spotAnti = 0.;
		
		Paths[i]->InsertValue(last_spot);
		PathsAntithetic[i]->InsertValue(last_spotAnti);
		
		for (size_t t = 0; t < nbSteps; ++t)
		{
			next_spot = last_spot * (1+ rate*dt+ sqrt(dt)*Brownian(i, t));
			next_spotAnti = last_spotAnti * (1 + rate*dt + sqrt(dt)*BrownianAntithetic(i, t));
			
			Paths[i]->InsertValue(next_spot);
			PathsAntithetic[i]->InsertValue(next_spotAnti);
			
			last_spot = next_spot;
			last_spotAnti = next_spotAnti;
		}
	}

};

