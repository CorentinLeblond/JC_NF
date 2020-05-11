#include "RandomNumbers.hpp"
#include <math.h>

#define MAXBIT 30
#define MAXDIM 6
//////////////////////////////////////////////////////////////////////////////////////////////
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
		// std::cout << "modula =  " << n % base << std::endl;
		// std::cout << "q =  " << q << std::endl;
		n /= base;
		// std::cout << "n = " << n << std::endl;
		bk /= base;
		// std::cout << "bk = " << bk << std::endl;
	}
	
	current+=1;
	return q;
};
////////////////////////////////////////////////////////////////////////////////
Sobol::Sobol(myLong inputseed): QuasiGenerator(inputseed)
{
	n = -1;
	sequence(n);
	n = 2;
	std::cout<<"Sobol Constructor"<<std::endl;
	requireNewSimulation = true;
	
};

double Sobol::generate()
{
	if(requireNewSimulation == true) 
	{
		// std::cout<<"Sobol Sequence generation"<<std::endl;
		sequence(n);
		// std::cout<<"Sobol Sequence generation Done"<<std::endl;
		requireNewSimulation = false;
		// std::cout<<"0D elemente: "<<x[0]<<std::endl;
		// std::cout<<"1D elemente: "<<x[1]<<std::endl;
		// std::cout<<"size: "<<x.size()<<std::endl;
		return x[1];
	}
	else
	{
		requireNewSimulation = true;
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
// int Sobol::GetN()
// {
	// return n;
// };
int IMIN(int a,int b)
{
	if(a <= b) return a;
	
	if(a > b) return b;
};
////////////////////////////////////////////////////////////////////////////////
// Kakutani::Kakutani(myLong inputseed):
	// QuasiGenerator(inputseed)
// {
// };
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
//Constructor, we use it to compute the covariance matrix thanks to inputs, it will be used in classes filles
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
	B = CovarMatrix.Cholesky();
	//std::cout << "Lower cholesky" <<std::endl;
	//LowerCholesly.Print();
};

matrix GaussianVectorCholesky::CorrelatedGaussianVector()
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

GaussianVectorDiag::GaussianVectorDiag(Normal* inputngnr, matrix inputSigma, matrix inputcorrel,matrix inputvarcovar)
	: GaussianVector(inputngnr, inputSigma,  inputcorrel, inputvarcovar)
{
	std::cout<<"GaussianVectorDiag Constructor"<<std::endl;
	v.Resize(CovarMatrix.nb_rows(),CovarMatrix.nb_rows());
	d.Resize(CovarMatrix.nb_rows(),1);
	n = CovarMatrix.nb_rows();
	nrot= 0;
	jacobi(CovarMatrix,n,d, v, nrot);
	std::cout<<"CovarMatrix after function called"<<std::endl;
	CovarMatrix.Print();
	std::cout<<"d matrix of eigevalues"<<std::endl;
	d.Print();
	std::cout<<"v matrix of eigenvectors"<<std::endl;
	v.Print();
	std::cout<<"nrot"<<nrot<<std::endl;
	d.Diagonalization();
	std::cout<<"d matrix of eigevalues"<<std::endl;
	d.Print();
	std::cout<<"d sqrt"<<std::endl;
	d.SQRT();
	d.Print();
	B = v*d;
	std::cout<<"output matrix"<<std::endl;
	B.Print();
};

#define ROTATE(a,i,j,k,l) g=a(i,j);h=a(k,l);a(i,j)=g-s*(h+g*tau);\
a(k,l)=h+s*(g-h*tau);

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
void jacobi(matrix& a, int n, matrix& d, matrix& v, int& nrot)
{/* 
	Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On
	output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
	v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
	a. nrot returns the number of Jacobi rotations that were required. */
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	b=vector(1,n);
	z=vector(1,n);
	for(ip=0;ip<n;ip++) 
	{
		for(iq=0;iq<n;iq++) v(ip,iq)=0.0;
		v(ip,ip)=1.0;
	}
	for (ip=0;ip<n;ip++) 
	{
		b[ip]=d(ip,0)=a(ip,ip);
		z[ip]=0.0;
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
			free_vector(z,1,n);
			free_vector(b,1,n);
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
					z[ip] -= h;
					z[iq] += h;
					d(ip,0) -= h;
					d(iq,0) += h;
					a(ip,iq)=0.0;
					for(j=0;j<ip-1;j++)
					{
						ROTATE(a,j,ip,j,iq)
					}
					for(j=ip;j<iq-1;j++)
					{
						ROTATE(a,ip,j,j,iq)
					}
					for(j=iq;j<n;j++)
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
			b[ip] += z[ip];
			d(ip,0)=b[ip];
			z[ip]=0.0; 
		}
	}
	std::cout << "Too many iterations in routine jacobi"<<std::endl;
}
#define NR_END 1
#define FREE_ARG char*
void free_vector(double *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
double *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	double *v;
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) std::cout<<"allocation failure in vector()"<<std::endl;
	return v-nl+NR_END;
}
////////////////////////////////////////////////////////////////////
// Python Version
		// ## module jacobi
	// ’’’ lam,x = jacobi(a,tol = 1.0e-9).
	// Solution of std. eigenvalue problem [a]{x} = lambda{x}
	// by Jacobi’s method. Returns eigenvalues in vector {lam}
	// and the eigenvectors as columns of matrix [x].
	// ’’’
	// from numarray import array,identity,diagonal
	// from math import sqrt
	// def jacobi(a,tol = 1.0e-9):
	// def maxElem(a): # Find largest off-diag. element a[k,l]
	// n = len(a)
	// aMax = 0.0
	
	// for i in range(n-1):
	// for j in range(i+1,n):
	// if abs(a[i,j]) >= aMax:
	// aMax = abs(a[i,j])
	// k = i; l = j
	// return aMax,k,l
	// def rotate(a,p,k,l): # Rotate to make a[k,l] = 0
	// n = len(a)
	// aDiff = a[l,l] - a[k,k]
	// if abs(a[k,l]) < abs(aDiff)*1.0e-36: t = a[k,l]/aDiff
	// else:
	// phi = aDiff/(2.0*a[k,l])
	// t = 1.0/(abs(phi) + sqrt(phi**2 + 1.0))
	// if phi < 0.0: t = -t
	// c = 1.0/sqrt(t**2 + 1.0); s = t*c
	// tau = s/(1.0 + c)
	// temp = a[k,l]
	// a[k,l] = 0.0
	// a[k,k] = a[k,k] - t*temp
	// a[l,l] = a[l,l] + t*temp
	// for i in range(k): # Case ofi<k
	// temp = a[i,k]
	// a[i,k] = temp - s*(a[i,l] + tau*temp)
	// a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
	// for i in range(k+1,l): # Case ofk<i<l
	// temp = a[k,i]
	// a[k,i] = temp - s*(a[i,l] + tau*a[k,i])
	// a[i,l] = a[i,l] + s*(temp - tau*a[i,l])
	// for i in range(l+1,n): # Case ofi>l
	// temp = a[k,i]
	// a[k,i] = temp - s*(a[l,i] + tau*temp)
	// a[l,i] = a[l,i] + s*(temp - tau*a[l,i])
	// for i in range(n): # Update transformation matrix
	// temp = p[i,k]
	// p[i,k] = temp - s*(p[i,l] + tau*p[i,k])
	// p[i,l] = p[i,l] + s*(temp - tau*p[i,l])
	// n = len(a)
	// maxRot = 5*(n**2)
	// p = identity(n)*1.0 # Initialize transformation matrix
	// for i in range(maxRot): # Jacobi rotation loop
	// aMax,k,l = maxElem(a)
	// if aMax < tol: return diagonal(a),p
	// rotate(a,p,k,l)
	// print ’Jacobi method did not converge’
// };

