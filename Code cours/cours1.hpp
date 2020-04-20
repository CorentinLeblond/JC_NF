
#include <iostream>
#include <vector>
//#include <cmath>
//#include <stdlib.h> 

//#include "continuous.hpp"

typedef unsigned unsigned long  myLong; 

class RandomGenerator
{
	public:
	
		RandomGenerator(){};
		~RandomGenerator(){};
		
		virtual double generate() = 0; //Pas de corps, les classes filles vont réimplémenter generate
		//Comme on ne redéfinit pas le calcul de moyenne pas besoin de mettre en virtual
		
		double mean(myLong nbsim);
		double var(myLong nbsim);
		
		myLong getCurrent(){return current;};
	
	protected:
	
		myLong current;
	
};

//Conteneur, pas d'attribut particulier
class UniformGenerator: public RandomGenerator
{
	public:
	
		UniformGenerator();
		~UniformGenerator(){};
		
		//double generate();
		//mean(unsigned long nbsim);

	
	
};
class PseudoGenerator: public UniformGenerator
{
	public:
	
		PseudoGenerator(myLong inputseed);
		PseudoGenerator(){};
		~PseudoGenerator(){};
		
		//double generate();
		//Bien de garder en mémoire le seed
		//mean(myLong nbsim);
	
	protected:
	
		myLong seed;
	
	
};
class LinearCongruential: public PseudoGenerator
{
	public:
	
	LinearCongruential(myLong inputmultiplier, myLong inputincrement, myLong inputmodulus, myLong inputseed);
	LinearCongruential(){};
	~LinearCongruential(){};
	
	double generate();
	//mean(myLong long nbsim);
	
	private:
	
	myLong multiplier;
	myLong increment;
	myLong modulus;
 
	
	
};
class EcuyerCombined: public PseudoGenerator
{
	public:
		
		EcuyerCombined();
		~EcuyerCombined(){};
		
		double generate();
		
	protected:
	
		LinearCongruential firstLinear;
		LinearCongruential secondLinear;
};


class DiscreteGenerator: public RandomGenerator
{
	public:
	
		DiscreteGenerator(){};
		~DiscreteGenerator(){};
	
	//double generate();
	// mean(unsigned long nbsim);
	
};

class HeadTail: public DiscreteGenerator
{
	public:
	
		HeadTail(UniformGenerator* inputugnr);
		~HeadTail(){};
		
		double generate();
	private:
	
		UniformGenerator* ugnr;
	
	
	//double generate();
	// mean(unsigned long nbsim);
};

class Bernoulli: public DiscreteGenerator
{
	public:
	
		Bernoulli(UniformGenerator* inputugnr,double inputp);
		~Bernoulli(){};
		
		double generate();
		
		
	
	private:
	
		UniformGenerator* ugnr;
		double p;
	
	
	//double generate();
	// mean(unsigned long nbsim);
};

class Binomial: public DiscreteGenerator
{
	public:
	
		Binomial(UniformGenerator* inputugnr,double inputp,myLong inputNbSimul);
		~Binomial(){};
		
		double generate();
		
		
	
	private:
	
		UniformGenerator* ugnr;
		double p;
		double nbSimul;
	
	
	//double generate();
	// mean(unsigned long nbsim);
};

class FiniteSet: public DiscreteGenerator
{
	public:
	
		FiniteSet(UniformGenerator* inputugnr,const std::vector<double>& inputproba);
		~FiniteSet(){};
		
		double generate();
	
	private:
	
		UniformGenerator* ugnr;
		std::vector<double> proba;
	
	
	//double generate();
	// mean(unsigned long nbsim);
};

class Poisson : public DiscreteGenerator
{
	protected:
	
		double Lambda;
		UniformGenerator* ugnr;
		

	public:
		Poisson(UniformGenerator* ugnr,double lambda);
		~Poisson(){};
		myLong Factorial(myLong n);
};

class PoissonFirstAlgo : public Poisson
{
	public:
		PoissonFirstAlgo(UniformGenerator* ugnr,double lambda);
		~PoissonFirstAlgo(){};
		
		double generate();
};

class ContinuousGenerator :public RandomGenerator
{
	public:
		ContinuousGenerator(){};
		~ContinuousGenerator(){};
};
class Exponential : public ContinuousGenerator
{
	protected:

		double Lambda;
		UniformGenerator* ugnr;

	public:
		Exponential(UniformGenerator* inputugnr,double inputLambda);
		~Exponential(){};
		
		double generate();
};

class Normal : public ContinuousGenerator
{
	protected:
		double Mu;
		double Sigma;
		UniformGenerator* ugnr;

	public:
		Normal(UniformGenerator* inputugnr,double inputMu, double inputSigma);
		~Normal();
};

class NormalBoxMuller : public Normal
{
	private:
		bool requireNewSimulation;
		double X;
		double Y;

	public:
		NormalBoxMuller(UniformGenerator* inputugnr,double inputMu, double inputSigma);
		double generate();
};

// class NormalCLT : public Normal
// {
	// public:
		// NormalCLT(double inputMu, double inputSigma);
		// double Generate();
// };

// class NormalRejectionSampling : public Normal
// {
	// public:
		// NormalRejectionSampling(double inputMu, double inputSigma);
		// double Generate();
// };
// class Poisson : public DiscreteGenerator
// {
	// protected:
		// double lambda;
		// double factorial();
		
	// public:
	
		// Poisson();
		// ~Poisson(){};
	
		
// };

// class PoissonFirstAlgo : public Poisson
// {
	// public:
	
		// PoissonFirstAlgo(double lambda);
		// double generate();
	
	
// };





	
	
	
