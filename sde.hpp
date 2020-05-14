#pragma once 
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "RandomNumbers.hpp"

/* This file contains all algorithms that are used to generate random paths for underlyings */


//Each single path object reprensents one path of one asset
class SinglePath
{
	public:
		
		SinglePath(double start, double end, size_t nbSteps);
		~SinglePath(){};
		
		void InsertValue(double val);
		double GetState(double time);
		std::vector<double> GetAllValues();
		
		
	protected:
	
		std::vector<double> Values;
		double StartTime;
		double EndTime;
		size_t Nbsteps;
	
};


class RandomProcess
{
	public:
		RandomProcess(){};
		RandomProcess(RandomGenerator* gen, int dim);
		~RandomProcess(){};
		
		virtual void Simulate(double startTime,double EndTime, size_t nbSteps) = 0;
		SinglePath* GetPath(int pos = 0);
		matrix GetAllPaths();
		matrix GetAllPathsAnti();
		double Expectation(size_t nb_simulation);
		double Get_Dt();
		double Get_rate();
		
	protected:
	
		RandomGenerator* Generator;
		std::vector<SinglePath*> Paths;
		int dimension;
		std::vector<SinglePath*> PathsAntithetic;
		matrix BrownianAntithetic;
		double dt_sde;
		double rate_sde;
		
};


class Brownian1D : public RandomProcess
{
	public:
	
		Brownian1D(RandomGenerator* Gen);
		~Brownian1D(){};
		
		void Simulate(double startTime,double endTime,size_t nbSteps);
};

class BlackScholes1D : public RandomProcess
{
	public:
	
		BlackScholes1D(RandomGenerator* Gen, double spot, double rate, double vol);
		~BlackScholes1D(){};
		
	protected:
	
		double Spot;
		double Rate;
		double Vol;
};

class BSEuler1D : public BlackScholes1D
{
	public:
	
		BSEuler1D(RandomGenerator* Gen, double spot, double rate, double vol);
		~BSEuler1D(){};
		
		void Simulate(double startTime,double EndTime,size_t nbSteps);
		
};

class BlackScholes2D:public RandomProcess
{
	public:
	
		BlackScholes2D(RandomGenerator* Gen, double spot1,double spot2,double rate1, double rate2
		,double vol1,double vol2,double rho);
		~BlackScholes2D(){};
		
	protected:
	
		double Spot1;
		double Rate1;
		double Vol1;
		double Spot2;
		double Rate2;
		double Vol2;
		double Rho;
};
class BSEuler2D : public BlackScholes2D
{
	public:
	
		BSEuler2D(RandomGenerator* Gen, double spot1,double spot2,double rate1, double rate2
		,double vol1,double vol2,double rho);
		~BSEuler2D(){};
		
		void Simulate(double startTime,double EndTime,size_t nbSteps);
		
};

class BlackScholesND : public RandomProcess
{
	public:

		BlackScholesND() {};
		BlackScholesND(GaussianVector* CorrelGaussian, matrix spot_vec,double inputrate);
		~BlackScholesND() {};

	protected:

		matrix V_spot;
		double rate;
		matrix Brownian;
		GaussianVector* m_gaussian;

};

class BSEulerND : public BlackScholesND

{
	public:
		BSEulerND() {};
		BSEulerND(GaussianVector* CorrelGaussian, matrix spot_vec,double inputrate);
		~BSEulerND() {};
		void Simulate(double startTime, double EndTime, size_t nbSteps);

};

class BSEulerNDAntithetic : public BlackScholesND
{
	public:
	
		BSEulerNDAntithetic() {};
		BSEulerNDAntithetic(GaussianVector* CorrelGaussian, matrix spot_vec,double inputrate);
		~BSEulerNDAntithetic() {};
		void Simulate(double startTime, double EndTime, size_t nbSteps);



};


