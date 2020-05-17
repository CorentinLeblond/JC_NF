#include "RandomNumbers.hpp"
#include <math.h>

#define MAXBIT 30
#define MAXDIM 6


/*This file contains all algorithms used to generate random distributions, from uniform to correlated gaussian vectors */

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
{};

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
{};

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
//////////////////////////////////////////////////////////////////////

/*Class used to generate quasi random monte carlo (Sobol and VanDerCorput)*/

QuasiGenerator::QuasiGenerator(myLong inputseed):
	seed(inputseed)
{
	current = inputseed;
};
VanDerCorput::VanDerCorput(int inputbase,myLong inputseed):
	QuasiGenerator(inputseed),
	base(inputbase)
{
};

double VanDerCorput::generate()
{
	n = current;
	q = 0.;
	bk = (double) 1/base;
	while(n>0)
	{
		q += (n % base)*bk;

		n /= base;

		bk /= base;

	}
	
	current+=1;
	return q;
};
////////////////////////////////////////////////////////////////////////////////
Sobol::Sobol(myLong inputseed): QuasiGenerator(inputseed)
{
	//We make use of the Sobol Sequence into the Box Muller Normal Generator
	//We generate a 2D Sobol sequence to obtain two independent distributions 
	
	//Call once with n as negative value to initialize the generation of the sequence
	n = -1;
	sequence(n);
	//Then n = 2 to set the 2D generation
	n = 2;
	std::cout<<"Sobol Constructor"<<std::endl;
	requireNewSimulation = true;
	
};

double Sobol::generate()
{
	//Each time the method is called, we generate a draw (2D) and return either the first dimension or the second one
	if(requireNewSimulation == true) 
	{
		sequence(n);

		requireNewSimulation = false;

		//Return 1D sequence
		return x[1];
	}
	else
	{
		requireNewSimulation = true;
		//Return 2nd Dimension result
		return x[2];
		
	}
};

void Sobol::sequence(int inputn)
{
	/*
	When n is negative, internally initializes a set of MAXBIT direction numbers for each of MAXDIM
	different Sobol’ sequences. When n is positive (but ≤MAXDIM), returns as the vector x[1..n]
	the next values from n of these sequences. (n must not be changed between initializations.)
	*/
// x.resize(n)
	int j,k,l;
	unsigned long i,im,ipp;
	static float fac;
	static unsigned long in,ix[MAXDIM+1],*iu[MAXBIT+1];
	static unsigned long mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
	static unsigned long ip[MAXDIM+1]={0,0,1,1,2,1,4};
	static unsigned long iv[MAXDIM*MAXBIT+1]={
	0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};
	if (n < 0) 
	{
		for (k=1;k<=MAXDIM;k++) ix[k]=0;
		in=0;
		if (iv[1] != 1) return;
		fac=1.0/(1L << MAXBIT);
		for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM) iu[j] = &iv[k];
		// To allowboth 1D and 2D addressing.
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
		x.resize(n+1);
		im=in++;
		for (j=1;j<=MAXBIT;j++) 
		{
			if (!(im & 1)) break;
			im >>= 1;
		}
		if (j > MAXBIT) std::cout<<("MAXBIT too small in sobseq")<<std::endl;
		im=(j-1)*MAXDIM;
		for (k=1;k<=IMIN(n,MAXDIM);k++)
		{
			ix[k] ^= iv[im+k];
			x[k]=ix[k]*fac;
		}
	}
};

int IMIN(int a,int b)
{
	if(a <= b) return a;
	
	if(a > b) return b;
};

//////////////////////////////////////////////////////////////////////
const long double PI = 3.141592653589793238L;

Normal::Normal(UniformGenerator* inputugnr,double inputMu, double inputSigma) 
	: Mu(inputMu),ugnr(inputugnr)
{
	if(inputSigma < 0)
		throw std::exception("The variance must be strictly positive for Normal distribution");
	Sigma = inputSigma;
}

Normal::~Normal()
{}

NormalBoxMuller::NormalBoxMuller(UniformGenerator* inputugnr,double inputMu, double inputSigma) 
	: 
	Normal(inputugnr,inputMu, inputSigma)
{
	requireNewSimulation = true;
}

double NormalBoxMuller::generate()
{
	if(requireNewSimulation == true)
	{
		double gnr1 = ugnr -> generate();
		double gnr2 = ugnr -> generate();
		
		
		// std::cout << "gnr1 " <<gnr1<<std::endl;
		// std::cout << "gnr2 " <<gnr2<<std::endl;
		
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
};
NormalBoxMullerVDC::NormalBoxMullerVDC(UniformGenerator* inputugnr,UniformGenerator* inputugnr2,double inputMu, double inputSigma) 
	: 
	NormalBoxMuller(inputugnr,inputMu, inputSigma),
	ugnr2(inputugnr2)
{
	requireNewSimulation = true;
}
double NormalBoxMullerVDC::generate()
{
	if(requireNewSimulation == true)
	{
		double gnr1 = ugnr -> generate();
		double gnr2 = ugnr2 -> generate();
		
		// std::cout << "gnr1 " <<gnr1<<std::endl;
		// std::cout << "gnr2 " <<gnr2<<std::endl;
		
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
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Correlated Gaussian Vector Generation

GaussianVector::GaussianVector(Normal* inputngnr, matrix inputSigma, matrix inputcorrel,matrix inputvarcovar):
	Ngnr(inputngnr),
	Sigma(inputSigma), //Vector of vol for all assets
	Correl_matrix(inputcorrel), //Matrix containing all correlation coefficients between -1,1
	CovarMatrix(inputvarcovar)
	{
		Nb_asset = inputSigma.nb_rows(); 
	};
	
matrix GaussianVector::GetB()
{
	return B;
};
GaussianVectorCholesky::GaussianVectorCholesky(Normal* inputngnr, matrix inputSigma, matrix inputcorrel,matrix inputvarcovar)
	: GaussianVector(inputngnr, inputSigma,  inputcorrel, inputvarcovar)
{
	//Generates the lower triangular matrix in the constructor from the method in matrix.cpp
	B = CovarMatrix.Cholesky();
};

matrix GaussianVectorCholesky::CorrelatedGaussianVector()
{
	matrix IndependentGaussian(Nb_asset,1);
	//We fill up our gaussian vector with independent draws 
	for(size_t i=0;i<Nb_asset;++i)
	{
		IndependentGaussian(i,0) = Ngnr->generate();
	}

	//Then, we only do the tranformation using the lowercholeskly and the Mu vector already available
	
	matrix a = B*IndependentGaussian;
	
	// Our independend Gaussian vector is now a correlated gaussian vector that we return
	
	return a;
};

GaussianVectorDiag::GaussianVectorDiag(Normal* inputngnr, matrix inputSigma, matrix inputcorrel,matrix inputvarcovar)
	: GaussianVector(inputngnr, inputSigma,  inputcorrel, inputvarcovar)
{
	
	// std::cout<<"GaussianVectorDiag Constructor"<<std::endl;
	v.Resize(CovarMatrix.nb_rows(),CovarMatrix.nb_rows());
	d.Resize(CovarMatrix.nb_rows(),1);
	n = CovarMatrix.nb_rows();
	nrot= 0;
	jacobi(CovarMatrix,n,d, v, nrot);
	// std::cout<<"CovarMatrix after function called"<<std::endl;
	// CovaqrMatrix.Print();
	// std::cout<<"d matrix of eigevalues"<<std::endl;
	// d.Print();
	// std::cout<<"v matrix of eigenvectors (normalized)"<<std::endl;
	// v.Print();
	// std::cout<<"nrotation "<<nrot<<std::endl;
	d.Diagonalization();
	d.SQRT();
	d.Print();
	B = v*d;
};



matrix GaussianVectorDiag::CorrelatedGaussianVector()
{
	matrix IndependentGaussian(Nb_asset,1);
	//We fill up our gaussian vector with independent draws 
	for(size_t i=0;i<Nb_asset;++i)
	{
		IndependentGaussian(i,0) = Ngnr->generate();
	}
	//std::cout << "IndependentGaussian1:"<<std::endl;
	//IndependentGaussian.Print();
	//std::cout << "IndependentGaussian ROWS: " << IndependentGaussian.nb_rows()<<std::endl;
	//std::cout << "IndependentGaussian COLS: " << IndependentGaussian.nb_cols()<<std::endl;
	//std::cout << "LowerC:"<<std::endl;

	//Then, we only do the tranformation using the lowercholeskly and the Mu vector already available
	
	matrix a = B*IndependentGaussian;
	
	//a.Print();
	// a += Nu;
	// Our independend Gaussian vector is now a correlated gaussian vector that we return
	
	return a;
};

#define ROTATE(a,i,j,k,l) g=a(i,j);h=a(k,l);a(i,j)=g-s*(h+g*tau);\
	a(k,l)=h+s*(g-h*tau);


//Diagonalization process for a symmetric matrix, we found that this method was the most efficient in our case where we work with the covariance matrix
//This method works when the matrix is symmetric and it is fairly fast
void jacobi(matrix& a, int n, matrix& d, matrix& v, int& nrot)
{/* 
	Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On
	output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
	v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
	a. nrot returns the number of Jacobi rotations that were required. */
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c;//,*b,*z;
	matrix b(n,1);
	matrix z(n,1);
	for(ip=0;ip<n;ip++) 
	{
		for(iq=0;iq<n;iq++) v(ip,iq)=0.0;
		v(ip,ip)=1.0;
	}
	for (ip=0;ip<n;ip++) 
	{
		b(ip,0)=d(ip,0)=a(ip,ip);
		z(ip,0)=0.0;
	}
	nrot=0;
	for (i=0;i<50;i++) 
	{
		sm=0.0;
		for(ip=0;ip<n-1;ip++)
		{
			for(iq=ip+1;iq<n;iq++)
				sm += fabs(a(ip,iq));
		}
		if (sm == 0.0) 
		{
			b.Clear();
			z.Clear();
			// free_vector(b,0,n);
			return;
		}
		if (i < 3)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for(ip=0;ip<n-1;ip++)
		{
			for(iq=ip+1;iq<n;iq++) 
			{
				g=100.0*fabs(a(ip,iq));
				if (i > 3 && (double)(fabs(d(ip,0))+g) == (double)fabs(d(ip,0))
					&& (double)(fabs(d(iq,0))+g) == (double)fabs(d(iq,0)))
					a(ip,iq)=0.0;
				else if (fabs(a(ip,iq)) > tresh)
				{
					h=d(iq,0)-d(ip,0);
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a(ip,iq))/h;
					else
					{
						theta=0.5*h/(a(ip,iq));
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a(ip,iq);
					z(ip,0) -= h;
					z(iq,0) += h;
					d(ip,0) -= h;
					d(iq,0) += h;
					a(ip,iq)=0.0;
					for(j=0;j<ip-1;j++)
					{
						ROTATE(a,j,ip,j,iq)
					}
					for(j=ip+1;j<iq-1;j++)
					{
						ROTATE(a,ip,j,j,iq)
					}
					for(j=iq+1;j<n;j++)
					{
						ROTATE(a,ip,j,iq,j)
					}
					for (j=0;j<n;j++)
					{
						ROTATE(v,j,ip,j,iq)
					}
					++nrot;
				}
			}
		}
		for (ip=0;ip<n;ip++)
		{
			b(ip,0) += z(ip,0);
			d(ip,0)=b(ip,0);
			z(ip,0)=0.0; 
		}
	}
};
