#include "cours1.hpp"
#include <math.h>

//RandomGenerator

double RandomGenerator::mean(myLong nbSim)
{

	double sum = 0;
	
	for(size_t i = 0; i < nbSim; ++i)
	{
		double unifn = generate();
		sum += unifn;
	}
	
	// std::cout << "mean via fonction mean : " << sum/nbSim << std::endl;
	return sum/nbSim;
};

double RandomGenerator::var(myLong nbSim)
{
	double sum = 0.;
	double _mean = 0.;
	
	_mean = mean(nbSim);
	
	double sum2 = 0.;
	
	for(size_t i = 0; i < nbSim; ++i)
	{
		double unifn = generate();
		sum2 += unifn*unifn;
	}
	
	sum2 /= nbSim;
	sum2 -= _mean*_mean;
	
	return sum2 ;
};

//UGenerator

UniformGenerator::UniformGenerator(){};

//PseudoGenerator

PseudoGenerator::PseudoGenerator(myLong inputseed):
	seed(inputseed)
{
	current = inputseed;
};

LinearCongruential::LinearCongruential(myLong inputmultiplier, myLong inputincrement, myLong inputmodulus,myLong inputseed):
	PseudoGenerator(inputseed),
	multiplier(inputmultiplier),
	increment(inputincrement),
	modulus(inputmodulus)
{
};

double LinearCongruential::generate()
{
	
	current = (multiplier*current + increment) % modulus;
	
	//Attention division entière, on fait un static cast en double pour obtenir un nb en tre 0 et 1
	return (double) current/modulus;
};

EcuyerCombined::EcuyerCombined():
	PseudoGenerator(),
	firstLinear(40014, 0, 2147483563,1),
	secondLinear(40692,0,2147483399,1)
{
};
double EcuyerCombined::generate()
{
	double unif1 = firstLinear.generate();
	double unif2 = secondLinear.generate();

	myLong x1 = firstLinear.getCurrent();
	myLong x2 = secondLinear.getCurrent();

	myLong m1 = 2147483563;

	current = (x1 - x2) % (m1 - 1);

	double result;

	if (current > 0)
	{
		result = (double) current / m1;
	}
	else if (current < 0)
	{
		result = (double) current / m1 + 1;
	}
	else
	{
		result = (double) (m1 - 1) / m1;
	}

	return result;
};

HeadTail::HeadTail(UniformGenerator* inputugnr):
	ugnr(inputugnr)	
{
};

double HeadTail::generate()
{
	double unifnb = ugnr->generate();
	
	if (unifnb>0.5)
		return 1;
	else
		return 0;
	
};

Bernoulli::Bernoulli(UniformGenerator* inputugnr, double inputp):
	ugnr(inputugnr),
	p(inputp)
{
};

double Bernoulli::generate()
{
	double unifnb = ugnr->generate();
	
	if (unifnb>p)
		return 1;
	else
		return 0;
	
};
Binomial::Binomial(UniformGenerator* inputugnr, double inputp,myLong inputn):
	ugnr(inputugnr),
	p(inputp),
	nbSimul(inputn)
{
};


double Binomial::generate()
{
	
	double compteur = 0.;
	
	for(size_t i = 0; i < nbSimul; ++i)
	{
		double unifnb = ugnr->generate();
		std::cout << "proba : " << unifnb << std::endl;
		
		if (unifnb>p)
			compteur+=1.;
		else
			compteur +=0.;
	}
	
	return compteur;
	
};

FiniteSet::FiniteSet(UniformGenerator* inputugnr,const std::vector<double>& inputproba):
	ugnr(inputugnr)
	// proba(inputproba)
	
{
	double sumProba = 0;
	for (size_t i = 0; i < inputproba.size(); ++i)
	{
		double currentProb = inputproba[i];
		if (currentProb < 0 || currentProb > 1)
		{
			throw std::exception("A probability must be between [0 , 1]");
		}
		else
		{
			sumProba += currentProb;
		}
	}
	if (sumProba != 1)
		
		throw std::exception("The sum of probabilities must be equal to 1");

	proba = inputproba;	
};

double FiniteSet::generate()
{
	size_t len_set = proba.size();
	std::cout << "len : " << len_set << std::endl;
	
	// double proba = 1./len_set;
	// std::cout << "p : " << proba << std::endl;
	
	double sumProba = 0;
	std::vector<double> probacumul{0.};
	
	for(size_t i = 0; i < len_set; ++i)
	{	
		double currentProb = proba[i];
		sumProba += currentProb;
		probacumul.push_back(sumProba);
	}	
	
	double unifnb = ugnr->generate();
	std::cout << "Tirage : " << unifnb << std::endl;

	for(size_t i = 0; i < len_set +1; ++i)
	{
		if ((unifnb>=probacumul[i]) && (unifnb < probacumul[i+1]))
		{
			std::cout << "Compris dans l'intervalle : " << probacumul[i]<< " et " <<probacumul[i+1]<< std::endl;
			return (double) (i+1);
			break; 
		}
	}
};

Poisson::Poisson(UniformGenerator* inputugnr,double lambda): 
	ugnr(inputugnr)
{
	if (lambda <= 0)
		throw std::exception("Lambda must be strictly positive for Poisson distribution");
	Lambda = lambda;
};
myLong Poisson::Factorial(myLong n)
{
	return (n == 1 || n == 0) ? 1 : Factorial(n - 1) * n;
};

PoissonFirstAlgo::PoissonFirstAlgo(UniformGenerator* ugnr,double lambda):
	Poisson(ugnr,lambda)
	{};
	
double PoissonFirstAlgo::generate()
{
	// double P = 0.;
	// double P_1 = 0.;
	
	myLong K = 0;
	
	double sum = 0.;

	do
	{
		K+=1;
		double gnr1 = ugnr -> generate();
		double expo =  -(log(gnr1)/Lambda);
		
		sum+=expo;
		
		// P_1 += exp(-Lambda)*pow(Lambda,K)/Factorial(K);
	}while(sum<=1);
	
	return (double) K-1;
	
};

Exponential::Exponential(UniformGenerator* inputugnr,double lambda):
	ugnr(inputugnr)
	{
		if (lambda <= 0)
			throw std::exception("Lambda must be strictly positive for Exponential distribution");
		
		Lambda = lambda;
	};
	
double Exponential::generate()
{
	double gnr1 = ugnr -> generate();
	
	return -(log(gnr1)/Lambda);
};

const long double PI = 3.141592653589793238L;

Normal::Normal(UniformGenerator* inputugnr,double inputMu, double inputSigma) : Mu(inputMu),ugnr(inputugnr)
{
	if(inputSigma < 0)
		throw std::exception("The variance must be strictly positive for Normal distribution");
	Sigma = inputSigma;
}

Normal::~Normal()
{}

NormalBoxMuller::NormalBoxMuller(UniformGenerator* inputugnr,double inputMu, double inputSigma) : Normal(inputugnr,inputMu, inputSigma)
{
	requireNewSimulation = true;
}

double NormalBoxMuller::generate()
{
	if(requireNewSimulation == true)
	{
		double gnr1 = ugnr -> generate();
		double gnr2 = ugnr -> generate();
		
		double R = sqrt(-2*log(gnr1));
		double O = 2*PI*gnr2;
		
		X = R*cos(O);
		Y = R*sin(O);
		
		X = X*Sigma + Mu;
		Y = Y*Sigma + Mu;
		
		requireNewSimulation = false;
		// std::cout << "X : " << X << std::endl;
		
		return X;
		
	}
	else
	{
		//std::cout << "Y : " << Y << std::endl;
		requireNewSimulation = true;
		return Y;

	}	
}
