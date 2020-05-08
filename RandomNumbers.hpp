
#include <iostream>
#include <vector>
#include <ostream>
#include "payoff.hpp"

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
//////////////////////////////////////////////////////////////////////
class QuasiGenerator: public UniformGenerator
{
	public:
	
		QuasiGenerator(myLong inputseed);
		QuasiGenerator(){};
		~QuasiGenerator(){};
		
		//double generate();
		//Bien de garder en mémoire le seed
		//mean(myLong nbsim);
	
	protected:
	
		myLong seed;
	
	
};
class VanDerCorput: public QuasiGenerator
{
	public:
	
		VanDerCorput(int inputbase = 2,myLong inputseed = 1);
		~VanDerCorput(){};
		
		double generate();
 
	private:
	
		int n;
		int base;
		double bk;
		double q;
	
};
class Sobol: public QuasiGenerator
{
	public:
	
		Sobol(myLong inputseed = 1);
		~Sobol(){};
		
		double generate();
		void sequence(int inputn);
		// int GetN();
		
	private:
		bool requireNewSimulation;
		std::vector<double> x;
		int n;
	
};
int IMIN(int a,int b);
// class Kakutani: public QuasiGenerator
// {
	// public:
	
		// Kakutani(myLong inputseed);
		// ~Kakutani(){};
		
		// double generate();
 
	// private:
	
		// int n;
		// int base;
		// double bk;
		// double q;
	
// };
//////////////////////////////////////////////////////////////////////

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
	protected:
		bool requireNewSimulation;
		double X;
		double Y;
	public:
		NormalBoxMuller(UniformGenerator* inputugnr,double inputMu, double inputSigma);
		~NormalBoxMuller(){};
		double generate();
};
class NormalBoxMullerVDC : public NormalBoxMuller
{
	protected:
		UniformGenerator* ugnr2;
	public:
		NormalBoxMullerVDC(UniformGenerator* inputugnr,UniformGenerator* inputugnr2,double inputMu, double inputSigma);
		~NormalBoxMullerVDC(){};
		double generate();
};

//On crée une classe mère "gaussian vector" pour permettre plusieurs techniques de décompositions via des classes filles
//Elle stocke les variables utiles qui seront utilisées dans les classes filles
//Deux potentielles classes filles pour deux types de décomposition (cholesly et l'autre méthode)
class GaussianVector
{
	public:
	
		GaussianVector(Normal* inputngnr,  matrix inputSigma, matrix inputcorrel,matrix varcovar);
		~GaussianVector(){};
		
		virtual matrix CorrelatedGaussianVector() =0;
		
	protected:
		
		size_t Nb_asset;
		matrix Sigma; //Vector of vol for all assets
		matrix Correl_matrix; //Matrix containing all correlation coefficients between -1,1
		matrix CovarMatrix;
		Normal* Ngnr;

};
class GaussianVectorCholesky: public GaussianVector
{
	public:
	
		GaussianVectorCholesky(Normal* inputngnr, matrix inputSigma, matrix inputcorrel,matrix varcovar);
		~GaussianVectorCholesky(){};
		
		// double generate();
				
		matrix CorrelatedGaussianVector();
		
	protected:
	
		matrix LowerCholesly;
};
class GaussianVectorDiag: public GaussianVector
{
	public:
	
		GaussianVectorDiag(Normal* inputngnr,  matrix inputSigma, matrix inputcorrel,matrix varcovar);
		~GaussianVectorDiag(){};
		
		// double generate();
		
		matrix CorrelatedGaussianVector();
		
	protected:
	
		matrix LowerDiag;
};


/*

int j,k,l;
unsigned long i,im,ipp;
static float fac;
static unsigned long in,ix[MAXDIM+1],*iu[MAXBIT+1];
static unsigned long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
static unsigned long ip[MAXDIM+1]={0,0,1,1,2,1,4};
static unsigned long iv[MAXDIM*MAXBIT+1]={
0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
if (*n < 0) 
{
	for (k=1;k<=MAXDIM;k++) ix[k]=0;
	in=0;
	if (iv[1] != 1) return;
	fac=1.0/(1L << MAXBIT);
	for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM) iu[j] = &iv[k];
	To allowboth 1D and 2D addressing.
	for (k=1;k<=MAXDIM;k++) 
	{
		for (j=1;j<=mdeg[k];j++) iu[j][k] <<= (MAXBIT-j);
		for (j=mdeg[k]+1;j<=MAXBIT;j++) 
		{
			ipp=ip[k];
			i=iu[j-mdeg[k]][k];
			i ^= (i >> mdeg[k]);
			for (l=mdeg[k]-1;l>=1;l--) 
			{
			if (ipp & 1) i ^= iu[j-l][k];
			ipp >>= 1;
			}
			iu[j][k]=i;
		}
	}
}
else 
{
	im=in++;
	for (j=1;j<=MAXBIT;j++) 
	{
		if (!(im & 1)) break;
		im >>= 1;
	}
	if (j > MAXBIT) nrerror("MAXBIT too small in sobseq");
	im=(j-1)*MAXDIM;
	for (k=1;k<=IMIN(*n,MAXDIM);k++)
	{
		ix[k] ^= iv[im+k];
		x[k]=ix[k]*fac;
	}
}
*/