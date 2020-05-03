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

matrix RandomProcess::GetAllPaths() {


	std::vector<double> element = Paths[0]->GetAllValues();
	matrix res(Paths.size(), element.size());
	for (size_t single = 0; single < Paths.size(); single++) {

		for (size_t i = 0; i < element.size(); i++) {

			res(single, i) = element[i];

		}

		element = Paths[single]->GetAllValues();

	}

	return res;

};

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


/////////////////////////////////////////////////////////////////////////////////////////////////////

BlackScholesND::BlackScholesND(Normal* N_gen, matrix spot_vec, matrix rate_vec, matrix Sigma_vec, matrix corr_matrix,
	matrix varcov)
	:m_Gen(N_gen), V_spot(spot_vec), V_Rate(rate_vec), V_vol(Sigma_vec), m_corr_matrix(corr_matrix), m_varcov(varcov),
	RandomProcess(m_Gen, V_spot.nb_rows())
{};

BSEulerND::BSEulerND(Normal* N_gen, matrix spot_vec, matrix rate_vec, matrix Sigma_vec, matrix corr_matrix,
	matrix varcov) :
	BlackScholesND(N_gen, spot_vec, rate_vec, Sigma_vec, corr_matrix, varcov)
{ };


void BSEulerND::Simulate(double startTime, double endTime, size_t nbSteps)
{
	Paths.clear();
	size_t assets = V_spot.nb_rows();
	double dt = (endTime - startTime) / nbSteps;
	Brownian.Resize(assets, nbSteps);
	matrix mean_vector(assets, 1);
	matrix ones(assets, 1);

	for (size_t i = 0; i < assets; ++i)
	{
		//create a vector full of one in order to apply a matrixwise 
		ones(i, 0) = 1.;
		//create the mean_vector at each time step of mu*dt
		mean_vector(i, 0) = V_Rate(i, 0) * dt;
	}

	for (size_t t = 0; t < nbSteps; ++t)
	{
		matrix X = GaussianVectorCholesky(m_Gen, mean_vector, V_vol, m_corr_matrix,
			m_varcov).CorrelatedGaussianVector();

		for (size_t i = 0; i < assets; ++i)
		{
			//create the matrix of all brownian motion 
			Brownian(i, t) = X(i, 0);
		}

	}

	SinglePath* Path = new SinglePath(startTime, endTime, nbSteps);
	for (size_t i = 0; i < assets; ++i)
	{
		double last_spot = V_spot(i, 0);
		double next_spot = 0.;
		Path->InsertValue(last_spot);
		for (size_t t = 0; t < nbSteps; ++t)
		{
			next_spot = last_spot * (1 + Brownian(i, t));
			Path->InsertValue(next_spot);
			last_spot = next_spot;

		}

		Paths.push_back(Path);
	}

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