#ifndef OMPRNG_H
#define OMPRNG_H

#include <omp.h>
#include "rngstream.h"
// #include "eigen/Eigen/Eigen"
#include "sys/time.h"
#include <cmath>

using namespace std;
// using namespace Eigen;


const double pi = std::acos(-1.0);



class omprng 
{
	private:
		int nprocs;
		RngStream *myRng;
		void randomSeed ();
	public:
		std::string name;
		omprng ();
		void fixedSeed (unsigned int);
		void setNumThreads (int);
		double runif ();
		double runif (double,double);
		double rnorm (double,double);
		double rexp (double);
		double rgamma (double,double);
		double rchisq (double);
		double rbeta (double,double);
		// void rmvn(const VectorXd,const MatrixXd,VectorXd& X);
};



omprng::omprng ()
{
	randomSeed();   				// default is random seed
	nprocs = omp_get_num_procs();	// determine number of processors
	myRng = new RngStream[nprocs];	// initialize RngStream
	myRng -> AdvanceState(3,3);		// advance state
}



void omprng::fixedSeed (unsigned int myInt) {
	unsigned long mySeed[6] = {myInt, myInt, myInt, myInt, myInt, myInt};
	RngStream::SetPackageSeed(mySeed);
	nprocs = omp_get_num_procs();
	myRng = new RngStream[nprocs];
	myRng -> AdvanceState(3,3);
}



void omprng::randomSeed () {
	timeval tim;
	gettimeofday(&tim, NULL);
	unsigned long int seed[6] = {(unsigned int)tim.tv_sec, (unsigned int)tim.tv_usec, (unsigned int) (tim.tv_sec+tim.tv_usec), 
		(unsigned int) abs(tim.tv_sec-tim.tv_usec), (unsigned int) abs(tim.tv_usec-tim.tv_sec), 5};
	RngStream::SetPackageSeed(seed);
}


void omprng::setNumThreads (int nt) {
	omp_set_num_threads(nt);		// set number of threads in OpenMP
	myRng = new RngStream[nt];		// initialize RngStream
}



double omprng::runif () 
{
	//generate random number between 0 and 1, i.e. from Unif(0,1)

	return(myRng[omp_get_thread_num()].RandU01());
}



double omprng::runif (double a, double b)
{
	//generate random number between a and b, i.e. from Unif(a,b)

	return(a + (b-a)*myRng[omp_get_thread_num()].RandU01());
}



double omprng::rnorm (double mu, double sigma) 
{
	// generate random number from X~N(mu,sigma)
	// f(x) = 1/sqrt(2*pi*si^2) * exp(-0.5*((x-mu)/si)^2)
	// -infty < x < infty, -infty < mu < infty, sigma > 0
	// E(X) = mu, Var(X) = sigma^2
	
	double u, v, z, x;

	u = runif();
	v = runif();

	z = sqrt(-2*log(u)) * cos(2*pi*v);
	x = mu + z*sigma;

	return(x);
}



double omprng::rexp (double theta) 
{
	// generate random number from X~exp(theta)
	// f(x) = 1 /theta * exp(-x/theta)
	// x > 0, theta > 0
	// E(X) = theta, Var(X) = theta^2

	double u, x;

	u = runif();
	x = -log(u)*theta;

	return(x);
}



double omprng::rgamma (double alpha, double beta) 
{
	// generate random number from X~gamma(alpha,beta)
	// f(x) = 1/(gamma(alpha)*beta^alpha) * x^(alpha-1) * exp(-x/beta)
	// x > 0, alpha > 0, beta > 0
	// E(X) = alpha*beta, Var(X) = alpha*beta^2

	double x=0.0;
	double delta, v0, v1, v2, v3, psi=1.0, nu=1.0e10;

	for(int i=0; i<floor(alpha); i++)
		x = x + rexp(1.0);

	if(alpha > floor(alpha)) { // alpha has fractional part
		delta = alpha - floor(alpha);
		while(nu > pow(psi,delta-1.0)*exp(-psi)) {
			v1 = runif(); v2 = runif(); v3 = runif();
			v0 = exp(1.0) / (exp(1.0) + delta);
			if(v1 <= v0) {
				psi = pow(v2, 1.0/delta);
				nu = v3 * pow(psi,delta-1.0);
			}
			else {
				psi = 1.0 - log(v2);
				nu = v3 * exp(-psi);		
			}
		}
		x = x + psi;
	}

	return(beta*x);
}



double omprng::rchisq (double df) 
{
	//generate random number from X~chisq(df)
	// f(x) = 1/(gamma(df/2)*2^(df/2)) * x^(df/2-1) * exp(-x/2)
	// x > 0, df = 1,2,3,...
	// E(X) = df, Var(X) = 2*df

	return(rgamma(df/2.0, 2.0));
}



double omprng::rbeta (double alpha, double beta)
{
	//generate random number from X~beta(alpha,beta)
	// f(x) = gamma(alpha+beta)/(gamma(alpha)*gamma(beta)) * x^(alpha-1) * (1-x)^(beta-1)
	// 0 < x < 1, alpha > 0, beta > 0
	// E(X) = alpha/(alpha+beta), Var(X) = alpha*beta / ((alpha+beta+1)*(alpha+beta)^2)

	double x, y;
	
	x = rgamma(alpha, 1.0);	
	y = rgamma(beta, 1.0);

	return(x/(x+y));
}



//The multivariate normal distribution is ALPHA code. Use with caution
// void omprng::rmvn(
//	const VectorXd mu, 
//	const MatrixXd var, 
//	VectorXd& X)
//{
//	int n = X.size();
//	VectorXd z(n);
//
//	// create object A to contain cholesky decomposition
//	Eigen::LLT<MatrixXd> A;
//
//	// generate z vector, where the z_i's are iid N(0,1)
//	for(int i=0; i<n; i++)
//		z(i) = rnorm(0.0,1.0);
//
//	// determine the mvn vector X = mu + Az
//	X = mu + A.compute(var).matrixL() * z;
//}


#endif
