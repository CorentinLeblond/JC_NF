#include "RandomNumbers.hpp"
#include <math.h>
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

// Kakutani::Kakutani(myLong inputseed):
	// QuasiGenerator(inputseed)
// {
// };
//////////////////////////////////////////////////////////////////////
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
};

//Constructor, we use it to compute the covariance matrix thanks to inputs, it will be used in classes filles
GaussianVector::GaussianVector(Normal* inputngnr,matrix inputMu, matrix inputSigma, matrix inputcorrel,matrix inputvarcovar):
	Ngnr(inputngnr),
	Nu(inputMu), //Vector of means for all assets
	Sigma(inputSigma), //Vector of vol for all assets
	Correl_matrix(inputcorrel), //Matrix containing all correlation coefficients between -1,1
	CovarMatrix(inputvarcovar)
	{
		Nb_asset = inputMu.nb_rows(); 
	};
	
	
GaussianVectorCholesky::GaussianVectorCholesky(Normal* inputngnr,matrix inputMu, matrix inputSigma, matrix inputcorrel,matrix inputvarcovar)
	: GaussianVector(inputngnr, inputMu,  inputSigma,  inputcorrel, inputvarcovar)
{
	LowerCholesly = CovarMatrix.Cholesky();
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
	
	matrix a = LowerCholesly*IndependentGaussian;
	
	//a.Print();
	a += Nu;
	// Our independend Gaussian vector is now a correlated gaussian vector that we return
	
	return a;
};

GaussianVectorDiag::GaussianVectorDiag(Normal* inputngnr,matrix inputMu, matrix inputSigma, matrix inputcorrel,matrix inputvarcovar)
	: GaussianVector(inputngnr, inputMu,  inputSigma,  inputcorrel, inputvarcovar)
{
////////////////////////////////////////
// C Version
	// void jacobi(float **a, int n, float d[], float **v, int *nrot)
	// Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On
	// output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
	// v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
	// a. nrot returns the number of Jacobi rotations that were required.
	// {
	// int j,iq,ip,i;
	// float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	// b=vector(1,n);
	// z=vector(1,n);
	// for (ip=1;ip<=n;ip++) { Initialize to the identity matrix.
	// for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
	// v[ip][ip]=1.0;
	// }
	// for (ip=1;ip<=n;ip++) { Initialize b and d to the diagonal
	// b[ip]=d[ip]=a[ip][ip]; of a.
	// z[ip]=0.0; This vector will accumulate terms
	// of the form tapq as in equation (11.1.14).
	// }
	// *nrot=0;
	// for (i=1;i<=50;i++) {
	// sm=0.0;
	// for (ip=1;ip<=n-1;ip++) { Sum off-diagonal elements.
	// for (iq=ip+1;iq<=n;iq++)
	// sm += fabs(a[ip][iq]);
	// }
	// if (sm == 0.0) { The normal return, which relies
	// on quadratic convergence to
	// machine underflow.
	// free_vector(z,1,n);
	// free_vector(b,1,n);
	// return;
	// }
	// if (i < 4)
	// tresh=0.2*sm/(n*n); ...on the first three sweeps.
	// else
	// tresh=0.0; ...thereafter.
	// for (ip=1;ip<=n-1;ip++) {
	// for (iq=ip+1;iq<=n;iq++) {
	// g=100.0*fabs(a[ip][iq]);
	// After four sweeps, skip the rotation if the off-diagonal element is small.
	// if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
	// && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
	// a[ip][iq]=0.0;
	// else if (fabs(a[ip][iq]) > tresh) {
	// h=d[iq]-d[ip];
	// if ((float)(fabs(h)+g) == (float)fabs(h))
	// t=(a[ip][iq])/h; t = 1/(2θ)
	// else {
	// theta=0.5*h/(a[ip][iq]); Equation (11.1.10).
	// t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	// if (theta < 0.0) t = -t;
	// }
	// c=1.0/sqrt(1+t*t);
	// s=t*c;
	// tau=s/(1.0+c);
	// h=t*a[ip][iq];
	// z[ip] -= h;
	// z[iq] += h;
	// d[ip] -= h;
	// d[iq] += h;
	// a[ip][iq]=0.0;
	// for (j=1;j<=ip-1;j++) { Case of rotations 1 ≤ j<p.
	// ROTATE(a,j,ip,j,iq)
	// }
	// for (j=ip+1;j<=iq-1;j++) { Case of rotations p<j<q.
	// ROTATE(a,ip,j,j,iq)
	// }
	// for (j=iq+1;j<=n;j++) { Case of rotations q<j ≤ n.
	// ROTATE(a,ip,j,iq,j)
	// }
	// for (j=1;j<=n;j++) {
	// ROTATE(v,j,ip,j,iq)
	// }
	// ++(*nrot);
	// }
	// }
	// }
	// for (ip=1;ip<=n;ip++) {
	// b[ip] += z[ip];
	// d[ip]=b[ip]; Update d with the sum of tapq,
	// z[ip]=0.0; and reinitialize z.
	// }
	// }
	// nrerror("Too many iterations in routine jacobi");
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
};

matrix GaussianVectorDiag::CorrelatedGaussianVector()
{
	matrix IndependentGaussian(Nb_asset,1);
	return IndependentGaussian;
};