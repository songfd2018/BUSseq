#ifdef _WIN32

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h> // pow, sqrt, lgamma
#include <cmath> //
#include "omprng.h"
# include <chrono> // time
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <direct.h>
#include <algorithm>    // sort
#include <vector>  
#include <R.h>
#include <Rinternals.h>

#endif

#ifdef linux

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h> // pow, sqrt, lgamma
#include <cmath> //
#include "omprng.h"
# include <chrono> // time
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <sys/stat.h>  //mkdir
#include <sys/types.h>
#include <algorithm>    // sort
#include <vector>  
#include <R.h>
#include <Rinternals.h>

#endif


using namespace std;

extern "C"{


void BUSseq_MCMC(int *Y_vec, int *Dim, int *seed, int *nc,
	int *iter_infor, char **dir_output, double *hyper, int *X_vec) ;

void BUSseq_inference(int *Y_vec, int *Dim,
	int *nc, int *iter_infor,
	char **dir_output,
	double *fdr_threshold,
	// posterior mean, mode and standard deviation
	double *alpha_est, double *alpha_sd,
	double * beta_est, double *beta_sd,
	double *nu_est, double *nu_sd,
	double *delta_est, double *delta_sd,
	double *gamma_est, double * gamma_sd,
	double *phi_est, double *phi_sd,
	double *pi_est, double *pi_sd,
	double *tau0_est, double *tau0_sd,
	double *p_est, double *p_sd,
	int *w_est, double *PPI_est,
	int *D_est, double *BIC);



int rand_cate(double* _prop, omprng _rng){
    int res = 0;
    double u = _rng.runif();
    while(u > _prop[res]){
        u = u - _prop[res];
        res++;
    }
    return res;
}

int rand_Ber(double _prob, omprng _rng){
    int res = 1;
    double u = _rng.runif();
    if(u > _prob){
        res = 0;
    }
    return res;
}

int rand_NB(double _r, double _mu, omprng _rng){
    
    double lambda, nu, p, u;
    double STEP = 500;// steps for large lambda in Poisson distribution
    int res = 0;

    // sample lambda from gamma distribution
    lambda = _rng.rgamma(_r, _mu/_r);
        
    // sample x from Poisson distribution
    nu = lambda;
    p = 1.0;

    do{
        res ++;
        u = _rng.runif();
        p = p * u;
        while(p < 1 & nu > 0){
            if(nu > STEP){
                p = p * exp(STEP);
                nu = nu - STEP;
            }else{
                p = p * exp(nu);
                nu = 0;
            }
        }
    }while(p > 1);

    res--;
    return res;
}

void rand_Dir(double *_xi, int _K, omprng _rng, double *_pi){

    double *rn_gam = new double[_K];
    double sum_gam = 0.0;
    for(int k = 0; k < _K; k++){
        rn_gam[k] = _rng.rgamma(_xi[k], 1.0);
        sum_gam += rn_gam[k];
    }
    for(int k = 0; k < _K; k++){
        _pi[k] = rn_gam[k] / sum_gam;
    }

    delete [] rn_gam;
}

void _update_logmu(int _B, int* _nb,
	int *_W, double _alpha, double* _beta, double* _nu, double* _delta,
	double* _logmu) {
	int cell_index = 0;
	int k;
	//auxiliry
	for (int b = 0; b < _B; b++) {
		for (int i = 0; i < _nb[b]; i++) {
			k = _W[cell_index];
			_logmu[cell_index] = _alpha + _beta[k] + _nu[b] + _delta[cell_index];
			cell_index++;
		}
	}
}

void _update_zx(int _B, int* _nb,
	double **_gamma, double *_phi, double *_logmu,
	int *_Y, omprng _rng,
	int *_X, int *_Z) {

	//int ind, ind_nu;
	//int ind_n = 0; //index the row of (b,i)	
	int cell_index = 0;
	double log_rat, u, temp_max, acc_rat;
	int temp_x;


	for (int b = 0; b < _B; b++) {
		for (int i = 0; i < _nb[b]; i++) {
			//ind = ind_n + j * _N;
			//ind_nu = b + j * _B;
			// printf("Sampling for cell %d.\n", cell_index + 1);
			// if(cell_index == 66){
			//     printf("Gamma_b0 = %f.\n", _gamma[b][0]);
			//     printf("mu_bi = %f\n", exp(_logmu[cell_index]));
			//     printf("phi = %f\n", _phi[b]);
			// }

			if (_Y[cell_index] == 0) {
				//update z_{big}
				if (_X[cell_index] == 0) {
					// printf("Sampling for Z.\n");
					log_rat = _gamma[b][0];
					_Z[cell_index] = rand_Ber(1 / (1 + exp(-log_rat)), _rng);
				}
				else {
					_Z[cell_index] = 1;
				}// end of if X == 0

				//Rprintf("z %d %d %d = %d\n", b, i, j, Dropout[index]);
				//update x_{big}
				if (_Z[cell_index] == 1) {// No dropout event
					// printf("Sampling for X.\n");
					//MH sampling
					temp_x = rand_NB(_phi[b], exp(_logmu[cell_index]), _rng);
					// printf("Sampling for X.\n");
					u = _rng.runif();
					//make the exponential be compuational
					temp_max = 0.0;
					if (temp_max < -_gamma[b][0] - _gamma[b][1] * temp_x) {
						temp_max = -_gamma[b][0] - _gamma[b][1] * temp_x;
					}
					if (temp_max < -_gamma[b][0] - _gamma[b][1] * _X[cell_index]) {
						temp_max = -_gamma[b][0] - _gamma[b][1] * _X[cell_index];
					}

					acc_rat = (exp(-temp_max) + exp(-_gamma[b][0] - _gamma[b][1] * _X[cell_index] - temp_max));
					acc_rat = acc_rat / (exp(-temp_max) + exp(-_gamma[b][0] - _gamma[b][1] * temp_x - temp_max));
					if (u < acc_rat) {
						_X[cell_index] = temp_x;
					}
				}
				else {
					_X[cell_index] = 0;
				}// end of if Z == 0
			}// end of if Y == 0
			cell_index++;
		}// end of i for
	}// end of b for
}

void _update_zx_optional(int _B, int* _nb, bool* drop,
	double **_gamma, double *_phi, double *_logmu,
	int *_Y, omprng _rng,
	int *_X, int *_Z) {

	//int ind, ind_nu;
	//int ind_n = 0; //index the row of (b,i)	
	int cell_index = 0;
	double log_rat, u, temp_max, acc_rat;
	int temp_x;


	for (int b = 0; b < _B; b++) {
		if (drop[b]) {
			for (int i = 0; i < _nb[b]; i++) {
				//ind = ind_n + j * _N;
				//ind_nu = b + j * _B;
				// printf("Sampling for cell %d.\n", cell_index + 1);
				// if(cell_index == 66){
				//     printf("Gamma_b0 = %f.\n", _gamma[b][0]);
				//     printf("mu_bi = %f\n", exp(_logmu[cell_index]));
				//     printf("phi = %f\n", _phi[b]);
				// }

				if (_Y[cell_index] == 0) {
					//update z_{big}
					if (_X[cell_index] == 0) {
						// printf("Sampling for Z.\n");
						log_rat = _gamma[b][0];
						_Z[cell_index] = rand_Ber(1 / (1 + exp(-log_rat)), _rng);
					}
					else {
						_Z[cell_index] = 1;
					}// end of if X == 0

					//Rprintf("z %d %d %d = %d\n", b, i, j, Dropout[index]);
					//update x_{big}
					if (_Z[cell_index] == 1) {// No dropout event
						// printf("Sampling for X.\n");
						//MH sampling
						temp_x = rand_NB(_phi[b], exp(_logmu[cell_index]), _rng);
						// printf("Sampling for X.\n");
						u = _rng.runif();
						//make the exponential be compuational
						temp_max = 0.0;
						if (temp_max < -_gamma[b][0] - _gamma[b][1] * temp_x) {
							temp_max = -_gamma[b][0] - _gamma[b][1] * temp_x;
						}
						if (temp_max < -_gamma[b][0] - _gamma[b][1] * _X[cell_index]) {
							temp_max = -_gamma[b][0] - _gamma[b][1] * _X[cell_index];
						}

						acc_rat = (exp(-temp_max) + exp(-_gamma[b][0] - _gamma[b][1] * _X[cell_index] - temp_max));
						acc_rat = acc_rat / (exp(-temp_max) + exp(-_gamma[b][0] - _gamma[b][1] * temp_x - temp_max));
						if (u < acc_rat) {
							_X[cell_index] = temp_x;
						}
					}
					else {
						_X[cell_index] = 0;
					}// end of if Z == 0
				}// end of if Y == 0
				cell_index++;
			}// end of i for
		}
		else {
			cell_index = cell_index + _nb[b];
		}
	}// end of b for
}

double _update_alpha(int _B, int* _nb,//dimension
	double _mu_a, double _sigma_a,//prior
	int *_W, double _alpha, double *_beta, double *_nu, double *_delta, double *_phi, //parameter
	int *_X, omprng _rng) {
	//int ind, ind_beta, ind_nu;
	//int ind_n; //index the row of (b,i)

	int cell_index = 0;

	//proposal
	double alpha_iter = _rng.rnorm(_alpha, 0.1);
	double logr = 0.0;
	int k;
	double res;

	//prior
	logr = logr - pow(alpha_iter - _mu_a, 2.0) / 2 / pow(_sigma_a, 2.0) + pow(_alpha - _mu_a, 2.0) / 2 / pow(_sigma_a, 2.0);

	for (int b = 0; b < _B; b++) {
		for (int i = 0; i < _nb[b]; i++) {
			k = _W[cell_index];

			//numerator
			logr = logr + alpha_iter * _X[cell_index] - (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(alpha_iter + _beta[k] + _nu[b] + _delta[cell_index]));
			//denomerator
			logr = logr - _alpha * _X[cell_index] + (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_alpha + _beta[k] + _nu[b] + _delta[cell_index]));

			cell_index++;
		}
	}

	if (logr > log(_rng.runif())) {
		res = alpha_iter;
	}
	else {
		res = _alpha;
	}

	return res;
}

void _update_l(int _K,//dimension
	double _p, double _tau0, double _tau1,//prior
	double *_beta, omprng _rng,//parameter
	int *_L) {
	//int ind_beta;
	double log_rat;
	for (int k = 1; k < _K; k++) {
		//ind_beta = j + k * _G;
		log_rat = 0.0; //the odds ratio of L_k = 1
		log_rat += log(_p) - log(1 - _p);
		log_rat += -log(_tau1) / 2.0 + log(_tau0) / 2.0;
		log_rat += -pow(_beta[k], 2) / 2.0 / _tau1;
		log_rat += pow(_beta[k], 2) / 2.0 / _tau0;
		_L[k] = rand_Ber(1.0 / (1.0 + exp(-log_rat)), _rng);
	}
}

void _update_beta(int _B, int *_nb, int _K,//dimension 
	double _tau0, double _tau1, int *_L,//prior
	int *_W, double _alpha, double *_nu, double *_delta, double *_phi, //parameter 	
	int *_X, omprng _rng,//latent variable
	double *_beta) {

	//index the row of (b,i)
	int cell_index, k;
	double *beta_iter = new double[_K];
	double *logr = new double[_K];

	for (k = 1; k < _K; k++) {
		//symmetric proposal
		beta_iter[k] = _rng.rnorm(_beta[k], 0.1);
		logr[k] = 0.0;

		//prior
		if (_L[k] == 1) {
			logr[k] = logr[k] - pow(beta_iter[k], 2.0) / 2 / _tau1 + pow(_beta[k], 2.0) / 2 / _tau1;
		}
		else {
			logr[k] = logr[k] - pow(beta_iter[k], 2.0) / 2 / _tau0 + pow(_beta[k], 2.0) / 2 / _tau0;
		}
	}

	cell_index = 0;
	for (int b = 0; b < _B; b++) {
		for (int i = 0; i < _nb[b]; i++) {
			k = _W[cell_index];
			//numerator
			logr[k] += beta_iter[k] * _X[cell_index];
			logr[k] += -(_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_alpha + beta_iter[k] + _nu[b] + _delta[cell_index]));
			//denomerator
			logr[k] += -_beta[k] * _X[cell_index];
			logr[k] += (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_alpha + _beta[k] + _nu[b] + _delta[cell_index]));

			cell_index++;
		}
	}

	for (k = 1; k < _K; k++) {
		if (logr[k] > log(_rng.runif())) {
			_beta[k] = beta_iter[k];
		}
	}
	delete[] logr;
	delete[] beta_iter;
}

void _update_nu(int _B, int *_nb,
	double *_mu_c, double _sigma_c,//prior
	int *_W, double _alpha, double *_beta, double *_delta, double *_phi,//parameter
	int *_X, omprng _rng, //latent variable
	double *_nu) {

	//int ind, ind_beta, ind_nu, ind_n;

	int cell_index = _nb[0];
	double nu_iter;
	double logr;
	int k;

	for (int b = 1; b < _B; b++) {
		//ind_nu = b + j * _B;

		//proposal
		nu_iter = _rng.rnorm(_nu[b], 0.1);
		logr = 0.0;

		//prior
		logr = logr - pow(nu_iter - _mu_c[b], 2.0) / 2 / pow(_sigma_c, 2.0) + pow(_nu[b] - _mu_c[b], 2.0) / 2 / pow(_sigma_c, 2.0);

		for (int i = 0; i < _nb[b]; i++) {

			//ind_beta = j + _w[ind_n] * _G;
			//ind = j * _N + ind_n;
			k = _W[cell_index];

			//numerator
			logr += nu_iter * _X[cell_index];
			logr += -(_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_alpha + _beta[k] + nu_iter + _delta[cell_index]));

			//denomerator
			logr += -_nu[b] * _X[cell_index];
			logr += (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_alpha + _beta[k] + _nu[b] + _delta[cell_index]));

			cell_index++;

			//ind_n = ind_n + 1;
		}

		if (logr > log(_rng.runif())) {
			_nu[b] = nu_iter;
		}
	}
}

void _update_phi(int _B, int *_nb,
	double *_phi_prior,//prior
	double *_logmu,//parameter
	int *_X, omprng _rng, //latent variable
	double *_phi) {
	//int ind, ind_nu, ind_n;
	int cell_index = 0;
	double phi_iter, logr;

	for (int b = 0; b < _B; b++) {

		phi_iter = _rng.rgamma(_phi[b], 1);
		logr = 0.0;

		for (int i = 0; i < _nb[b]; i++) {

			//numerator
			logr += lgamma(phi_iter + _X[cell_index]) - lgamma(phi_iter);
			logr += phi_iter * log(phi_iter) - (phi_iter + _X[cell_index]) * log(phi_iter + exp(_logmu[cell_index]));

			//denomerator
			logr += -lgamma(_phi[b] + _X[cell_index]) + lgamma(_phi[b]);
			logr += -_phi[b] * log(_phi[b]) + (_phi[b] + _X[cell_index]) * log(_phi[b] + exp(_logmu[cell_index]));

			cell_index++;
		}

		//prior
		logr += (_phi_prior[0] - 1) * log(phi_iter) - _phi_prior[1] * phi_iter;
		logr += -(_phi_prior[0] - 1) * log(_phi[b]) + _phi_prior[1] * _phi[b];

		//proposal
		logr += (phi_iter - 1.0) * log(_phi[b]) + phi_iter - lgamma(phi_iter);
		logr += -(_phi[b] - 1.0) * log(phi_iter) - _phi[b] + lgamma(_phi[b]);

		if (logr > log(_rng.runif())) {
			_phi[b] = phi_iter;
		}
	}
}

//Max value in a vector
double vec_max(double* vec,int n){
	double res = vec[0];
	for (int i = 1; i < n; i++) {
		if (res < vec[i]) {
		    res = vec[i];
		}
	}
	return res;
}

// FDR calculation
double fdrDEindicator(double *_PPI, double _kappa, int _G, int _K){
    
    double xi, fdr;
    double sum_xi = 0.0;
    int count_intri = 0;
    if(_kappa > 0){
        for(int g = 0; g < _G; g++){
            for(int k = 1; k < _K; k++){
                xi = 1 - _PPI[g * _K + k];
                if(xi <= _kappa){
                    sum_xi += xi;
                    count_intri ++;
                }
            }
        }
        fdr = sum_xi/count_intri;
    }else{
        fdr = 0.0;
    }
        
    return(fdr);
}

bool descreasing(double i,double j) { return (i>j); }

// calculate the DE posterior probability threshold
double postprob_DE_thr_fun(double *_PPI, double _fdr_threshold, int _G, int _K){
    
    double kappa;
    double postprob_DE_thr = 0.5;
    double fdr = 0.0;
    vector<double> vec_PPI;
    for(int g = 0; g < _G; g++){
        copy(_PPI + g * _K +1, _PPI + g * _K  + _K, back_inserter(vec_PPI));
    }
    // sort all PPIs decreasingly
    sort(vec_PPI.begin(),vec_PPI.end(),descreasing);
	 
    // unique values in the PPI vector
    vector<double>::iterator it;
    it = unique (vec_PPI.begin(), vec_PPI.end());                                                           //                ^
    vec_PPI.resize(distance(vec_PPI.begin(),it));

    int i = 0;
	kappa = 1 - vec_PPI[0];
	fdr = fdrDEindicator(_PPI, kappa, _G, _K);

    while(fdr <= _fdr_threshold & kappa <= postprob_DE_thr &i < (unsigned)vec_PPI.size()){
        i++;
		kappa = 1 - vec_PPI[i];
        fdr = fdrDEindicator(_PPI, kappa, _G, _K);
    }

	if(i == 0){
		kappa = 0.0;
	}else{
		kappa = 1 - vec_PPI[i-1];
	}
    
    return kappa;
}

int IG_index(double *_PPI, double _postprob_DE_thr, int _K){

    int IG_indicator = 0;
    for(int k = 1; k < _K; k++){
        if(_PPI[k] >= _postprob_DE_thr){
            IG_indicator = 1;
        }
    }
    
    return IG_indicator;
}

void BUSseq_MCMC(int *Y_vec, int *Dim, int *seed, int *nc,
	int *iter_infor, char **dir_output, double *hyper, int *X_vec) {
	////////////////////////////////////////
	//  0. Deal with the input from R  //
	////////////////////////////////////////
	// Load the dimension information and raw count data
	int N = Dim[0];
	int G = Dim[1];
	int B = Dim[2];
	int K = Dim[3];
	int *nb = &(Dim[4]);
	int *di = &(Dim[4 + B]);
	bool *Drop_ind = new bool[B];
	bool All_Drop = true;
	for (int b = 0; b < B; b++) {
		if (di[b] == 1) {
			Drop_ind[b] = true;
		}
		else if (di[b] == 0) {
			Drop_ind[b] = false;
			All_Drop = false;
		}
	}
	int **Y = new int*[G];
	for (int g = 0; g < G; g++) {
		Y[g] = &Y_vec[g * N];// for parallel G
	}
	// Set the seed for RNG
	omprng MCMC_Rng;
	MCMC_Rng.fixedSeed(seed[0]);
	// Set the number of cores for parallel
	omp_set_num_threads(nc[0]);
	// Set the number of iterations
	int iter_max = iter_infor[0]; // the overall iteration number 
	int iter_out = iter_infor[1]; // the number of iterations to print the posterior sampling 
	int iter_noupdate = iter_infor[2];
	// Set the directory to output the posterior sampling
	string output_dir(dir_output[0]);
	output_dir = output_dir + "/";

#ifdef linux
int check = mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

#ifdef _WIN32
if (0 != access(output_dir.c_str(), 0))
    {
        // if this folder not exist, create a new one.
        mkdir(output_dir.c_str()); 
    }
#endif
	
	string out_file;
	// file variable to output the posterior sampling
	ofstream post_File;
	// hyper-parameters
	// for prop Dir prior
	double xi = hyper[0];
	// for gamma_b0 Normal prior
	double sigma_z = hyper[1];
	// for gamma_b1 Gamma prior
	double *gamma_prior = &(hyper[2]);
	// for alpha Normal prior
	double sigma_a = hyper[4];
	// for beta spike and slab prior
	double tau1 = hyper[5];
	// for p beta prior
	double *p_prior = &(hyper[6]);
	// for tau0sq inverse-gamma prior
	double *tau0_prior = &(hyper[8]);
	// for nu Normal prior
	double sigma_c = hyper[10];
	// for delta Normal prior
	double sigma_d = hyper[11];
	// for phi Gamma prior
	double *phi_prior = &(hyper[12]);
	auto start_overall = chrono::system_clock::now();
	////////////////////////////////
	// 1. Initialize parameters //
	////////////////////////////////
	// cout << "Set initial values." << endl;
	// allocate memory for parameters
	// cell type proportion
	double temp_sum = (K + 1) * K / 2;// sum up 1 to K
	double **prop = new double*[B];
	for (int b = 0; b < B; b++) {
		prop[b] = new double[K];
		for (int k = 0; k < K; k++) {
			prop[b][k] = (K - k) / temp_sum;
		}
	}
	// cell type indicator
	int *W = new int[N];
	int *pointer_w = W;
	for (int b = 0; b < B; b++) {
		for (int i = 0; i < nb[b]; i++) {
			*pointer_w = rand_cate(prop[b], MCMC_Rng);
			pointer_w++;
		}
	}
	// cell-specific effect
	double *delta = new double[N];
	double *mu_d = new double[N];
	double *pointer_del = delta;
	int cell_index = 0;
	for (int b = 0; b < B; b++) {
		double sum_y0;
		for (int i = 0; i < nb[b]; i++) {
			if (i == 0) {
				*pointer_del = 0.0;// the first cell in each batch as reference
				sum_y0 = 0.0;// sum_g log(1+Y_b0g)
				for (int g = 0; g < G; g++) {
					sum_y0 += Y[g][cell_index];
				}
				pointer_del++;
			}
			else {
				double sum_y;
				sum_y = 0.0;// sum_g log(1+Y_big)
				for (int g = 0; g < G; g++) {
					sum_y += Y[g][cell_index];
				}
				*pointer_del = log(sum_y) - log(sum_y0);
				mu_d[cell_index] = *pointer_del;
				pointer_del++;
			}
			cell_index++;
		}
	}
	// baseline effects
	double *alpha = new double[G];
	double *mu_a = new double[G];
	// cell-type specific effects
	double **beta = new double*[G];
	for (int g = 0; g < G; g++) {
		beta[g] = new double[K];
	}
	// cout << "Calculate the empirical mean" << endl;
	for (int g = 0; g < G; g++) {
		double *sum_logy = new double[K];
		int *count_type = new int[K];
		for (int k = 0; k < K; k++) {
			sum_logy[k] = 0.0;
			count_type[k] = 0;
		}
		cell_index = 0;
		// only in the first batch
		for (int i = 0; i < nb[0]; i++) {
			int k = W[cell_index];
			sum_logy[k] += log(1 + Y[g][cell_index] / exp(delta[cell_index]));
			count_type[k] ++;
			cell_index++;
		}
		alpha[g] = sum_logy[0] / (count_type[0] + 1);// + 1 to prevent count_type[0] = 1
		mu_a[g] = alpha[g];
		beta[g][0] = 0.0;
		for (int k = 1; k < K; k++) {
			beta[g][k] = sum_logy[k] / (count_type[k] + 1) - alpha[g];
		}
		delete[] sum_logy;
		delete[] count_type;
	}
	// batch effects
	double **nu = new double*[G];
	double *mu_c = new double[B];
	for (int b = 0; b < B; b++) {
		mu_c[b] = 0.0;
	}
	for (int g = 0; g < G; g++) {
		nu[g] = new double[B];
		nu[g][0] = 0.0;
		double sum_logy0;
		sum_logy0 = 0.0;
		cell_index = 0;
		for (int i = 0; i < nb[0]; i++) {
			sum_logy0 += log(1 + Y[g][cell_index] / exp(delta[cell_index]));
			cell_index++;
		}
		for (int b = 1; b < B; b++) {
			double sum_logy = 0.0;
			for (int i = 0; i < nb[b]; i++) {
				sum_logy += log(1 + Y[g][cell_index] / exp(delta[cell_index]));
				cell_index++;
			}
			nu[g][b] = sum_logy / nb[b] - sum_logy0 / nb[0];
			mu_c[b] += nu[g][b];
		}

	}
	for (int b = 0; b < B; b++) {
		mu_c[b] = mu_c[b] / G;
	}
	// over-dispersion parameters
	double **phi = new double*[G];
	for (int g = 0; g < G; g++) {
		phi[g] = new double[B];
		for (int b = 0; b < B; b++) {
			phi[g][b] = 5.0;
		}
	}
	// intercept and odds ratio of dropout events
	double **gamma = new double*[B];
	for (int b = 0; b < B; b++) {
		gamma[b] = new double[2];
		gamma[b][0] = 0.0;
		if (Drop_ind[b]) {
			gamma[b][1] = -0.1;
		}
		else {
			gamma[b][1] = 0.0;
		}
	}
	// underlying true read count
	int **X = new int*[G];
	for (int g = 0; g < G; g++) {
		X[g] = &X_vec[g * N];
		for (int i = 0; i < N; i++) {
			X[g][i] = Y[g][i];
		}
	}
	// dropout indicator
	int **Z = new int*[G];
	for (int g = 0; g < G; g++) {
		Z[g] = new int[N];
		for (int i = 0; i < N; i++) {
			Z[g][i] = 0;
		}
	}
	// intrinsic gene proportion and variance of non-intrinsic genes
	double p, tau0; // p and tau0 as well as tau1
	p = MCMC_Rng.runif(0, 0.5);
	tau0 = 0.005;
	// intrinsic gene indicators
	int **L = new int*[G];
#pragma omp parallel for
	for (int g = 0; g < G; g++) {
		L[g] = new int[K];
		L[g][0] = 0;
		for (int k = 1; k < K; k++) {
			double log_rat;  // log odds of gene g in the cell type k being an intrinsic gene
			log_rat = log(p / (1 - p));
			log_rat = log_rat - log(tau1) / 2 - pow(beta[g][k], 2) / 2 / tau1;
			log_rat = log_rat + log(tau0) / 2 + pow(beta[g][k], 2) / 2 / tau0;
			L[g][k] = rand_Ber(1 / (1 + exp(-log_rat)), MCMC_Rng);
		}
	}
	// cout << "Allocate memory for recording the posterior sampling." << endl;
	// Allocate memory to store posterior sampling
	double **alpha_post = new double*[iter_out];
	for (int iter = 0; iter < iter_out; iter++) {
		alpha_post[iter] = new double[G];
	}
	double **beta_post = new double*[iter_out];
	for (int iter = 0; iter < iter_out; iter++) {
		beta_post[iter] = new double[G * K];
	}
	double **nu_post = new double*[iter_out];
	for (int iter = 0; iter < iter_out; iter++) {
		nu_post[iter] = new double[G * B];
	}
	double **delta_post = new double*[iter_out];
	for (int iter = 0; iter < iter_out; iter++) {
		delta_post[iter] = new double[N];
	}
	double **gamma_post = new double*[iter_out];
	for (int iter = 0; iter < iter_out; iter++) {
		gamma_post[iter] = new double[B * 2];
	}
	double **phi_post = new double*[iter_out];
	for (int iter = 0; iter < iter_out; iter++) {
		phi_post[iter] = new double[G * B];
	}
	double **pi_post = new double*[iter_out];
	for (int iter = 0; iter < iter_out; iter++) {
		pi_post[iter] = new double[B * K];
	}
	int **w_post = new int*[iter_out];
	for (int iter = 0; iter < iter_out; iter++) {
		w_post[iter] = new int[N];
	}
	double *p_post = new double[iter_out];
	double *tau0_post = new double[iter_out];
	int **l_post = new int*[iter_out];
	for (int iter = 0; iter < iter_out; iter++) {
		l_post[iter] = new int[G * K];
	}
	///////////////////////////
	// 2. MCMC sampling //
	///////////////////////////
	// auxiliary variable
	// cout << "Start MCMC sampling." << endl;
	double **logmu = new double*[G];
	for (int g = 0; g < G; g++) {
		logmu[g] = new double[N];
	}
#pragma omp parallel for
	for (int g = 0; g < G; g++) {
		_update_logmu(B, nb,
			W, alpha[g], beta[g], nu[g], delta,//parameter
			logmu[g]);
	}
	// divide the overall iterations into several parts with length iter_out
	// to control the RAM occupation
	int out_times = (iter_max - 1) / iter_out + 1;
	int iter_last = iter_max - iter_out * (out_times - 1);
	// cout << "Output the posterior sampling to txt file per " << iter_out << " iterations." << endl;
	int IND_UPDATE_PTAU0 = 0;
	// auxiliary variable to update gamma
	double gamma_iter, logr;//temp_pres, temp_iter, logr;
	double sigma_zsq = pow(sigma_z, 2.0);
	// auxiliary variable to update p and tau0
	double *p_postdist = new double[2];
	double *tau0_postdist = new double[2];
	// auxiliary variable  to update delta
	double delta_iter;
	// auxiliary variable to update w
	double *proposal_pi = new double[K];
	for (int k = 0; k < K; k++) {
		proposal_pi[k] = 1.0 / K;
	}
	int w_proposal, w_current;
	double log_proposal, log_current;
	// auxiliary variable to update prop
	double *count_w = new double[K];
	// auxiliary variable to store posterior sampling
	int q;

 	// variables for progress bar
	double progress = 0.0;
	double tot_iter = 0.0;
	int barWidth = 40;

	for (int t = 0; t < out_times; t++) {
		int iter_num = iter_out;
		if (t == out_times - 1) {
			iter_num = iter_last;
		}
		if (iter_noupdate > iter_num) {
			iter_noupdate = iter_noupdate - iter_num;
		}
		// cout << "Starting " << t * iter_out + 1 << "-" << t * iter_out + iter_num << " iterations for the " << t + 1 << "-th output." << endl;
		auto start_MCMC = chrono::system_clock::now();
		for (int iter = 0; iter < iter_num; iter++) {
			if (All_Drop) {
				/////////////////////////////////////
				//  1) update z_{big} and x_{big}  //
#pragma omp parallel for
				for (int g = 0; g < G; g++) {
					_update_zx(B, nb,
						gamma, phi[g], logmu[g],
						Y[g], MCMC_Rng, X[g], Z[g]);
				}
				/////////////////////////////////////////////////////
				//  2) update gamma_{b0} and gamma_{b1} // 

				// gamma_{b0}        
				cell_index = 0;
				for (int b = 0; b < B; b++) {
					gamma_iter = MCMC_Rng.rnorm(gamma[b][0], 0.1);

					//prior
					logr = -pow(gamma_iter, 2.0) / 2 / sigma_zsq + pow(gamma[b][0], 2.0) / 2 / sigma_zsq;

					for (int i = 0; i < nb[b]; i++) {
#pragma omp parallel
						{
							// auxiliary variable to update gamma
							double temp_pres, temp_iter;
							double logr_thread = 0.0;
							int* X_thread, *Z_thread;
#pragma omp for
							for (int g = 0; g < G; g++) {
								X_thread = X[g];
								Z_thread = Z[g];
								//numerator
								temp_iter = gamma_iter + X_thread[cell_index] * gamma[b][1];//prevent temp_iter is extremely large
								if (temp_iter > 0) {
									logr_thread += gamma_iter * Z_thread[cell_index] - temp_iter - log(1 + exp(-temp_iter));
								}
								else {
									logr_thread += gamma_iter * Z_thread[cell_index] - log(1 + exp(temp_iter));
								}

								//denomerator
								temp_pres = gamma[b][0] + X_thread[cell_index] * gamma[b][1];
								if (temp_pres > 0) {
									logr_thread += -gamma[b][0] * Z_thread[cell_index] + temp_pres + log(1 + exp(-temp_pres));
								}
								else {
									logr_thread += -gamma[b][0] * Z_thread[cell_index] + log(1 + exp(temp_pres));
								}
							}
#pragma omp critical
							{
								logr += logr_thread;
							}
						}// end of omp parallel
						cell_index++;
					}// end of i

					// cout << " logr = " << logr << endl;
					if (logr > log(MCMC_Rng.runif())) {
						gamma[b][0] = gamma_iter;
					}
				}// end of b

				// gamma_{b1}
				cell_index = 0;
				for (int b = 0; b < B; b++) {
					//pro posal
					gamma_iter = -MCMC_Rng.rgamma(-10 * gamma[b][1], 0.1);
					//prior
					logr = (gamma_prior[0] - 1) * (log(gamma_iter / gamma[b][1])) + gamma_prior[1] * (gamma_iter - gamma[b][1]);
					//proposal
					//numerator
					logr = logr - lgamma(-10 * gamma_iter) + (-10 * gamma_iter - 1) * log(-gamma[b][1]) - 10 * gamma_iter * log(10) + 10 * gamma[b][1];
					//denomerator
					logr = logr + lgamma(-10 * gamma[b][1]) - (-10 * gamma[b][1] - 1) * log(-gamma_iter) + 10 * gamma[b][1] * log(10) - 10 * gamma_iter;

					for (int i = 0; i < nb[b]; i++) {

#pragma omp parallel
						{
							// auxiliary variable to update gamma
							double temp_pres, temp_iter;
							double logr_thread = 0.0;
							int* X_thread, *Z_thread;
#pragma omp for
							for (int g = 0; g < G; g++) {
								X_thread = X[g];
								Z_thread = Z[g];
								//numerator
								temp_iter = gamma[b][0] + X_thread[cell_index] * gamma_iter;//prevent temp_iter is extremely large
								if (temp_iter > 0) {
									logr_thread += X_thread[cell_index] * gamma_iter * Z_thread[cell_index] - temp_iter - log(1 + exp(-temp_iter));
								}
								else {
									logr_thread += X_thread[cell_index] * gamma_iter * Z_thread[cell_index] - log(1 + exp(temp_iter));
								}
								//denomerator
								temp_pres = gamma[b][0] + X_thread[cell_index] * gamma[b][1];
								if (temp_pres > 0) {
									logr_thread += -X_thread[cell_index] * gamma[b][1] * Z_thread[cell_index] + temp_pres + log(1 + exp(-temp_pres));
								}
								else {
									logr_thread += -X_thread[cell_index] * gamma[b][1] * Z_thread[cell_index] + log(1 + exp(temp_pres));
								}
							}
#pragma omp critical
							{
								logr += logr_thread;
							}
						}
						cell_index++;
					}//end of i
					if (logr > log(MCMC_Rng.runif())) {
						gamma[b][1] = gamma_iter;
					}
				}// end of b
			}
			else {
				/////////////////////////////////////
				//  1) update z_{big} and x_{big}  //
				// auto start_zx = chrono::system_clock::now();
#pragma omp parallel for
				for (int g = 0; g < G; g++) {
					_update_zx_optional(B, nb, Drop_ind,
						gamma, phi[g], logmu[g],
						Y[g], MCMC_Rng,
						X[g], Z[g]);
				}
				//////////////////////////////////////////////////////
				//  2) update gamma_{b0} and gamma_{b1} // 
				// gamma_{b0}        
				cell_index = 0;
				for (int b = 0; b < B; b++) {
					if (Drop_ind[b]) {
						gamma_iter = MCMC_Rng.rnorm(gamma[b][0], 0.1);
						//prior
						logr = -pow(gamma_iter, 2.0) / 2 / sigma_zsq + pow(gamma[b][0], 2.0) / 2 / sigma_zsq;
						for (int i = 0; i < nb[b]; i++) {
#pragma omp parallel
							{
								// auxiliary variable to update gamma
								double temp_pres, temp_iter;
								double logr_thread = 0.0;
								int* X_thread, *Z_thread;
#pragma omp for
								for (int g = 0; g < G; g++) {
									X_thread = X[g];
									Z_thread = Z[g];
									//numerator
									temp_iter = gamma_iter + X_thread[cell_index] * gamma[b][1];//prevent temp_iter is extremely large
									if (temp_iter > 0) {
										logr_thread += gamma_iter * Z_thread[cell_index] - temp_iter - log(1 + exp(-temp_iter));
									}
									else {
										logr_thread += gamma_iter * Z_thread[cell_index] - log(1 + exp(temp_iter));
									}
									//denomerator
									temp_pres = gamma[b][0] + X_thread[cell_index] * gamma[b][1];
									if (temp_pres > 0) {
										logr_thread += -gamma[b][0] * Z_thread[cell_index] + temp_pres + log(1 + exp(-temp_pres));
									}
									else {
										logr_thread += -gamma[b][0] * Z_thread[cell_index] + log(1 + exp(temp_pres));
									}
								}
#pragma omp critical
								{
									logr += logr_thread;
								}
							}// end of omp parallel
							cell_index++;
						}// end of i
												// cout << " logr = " << logr << endl;
						if (logr > log(MCMC_Rng.runif())) {
							gamma[b][0] = gamma_iter;
						}
						else {
							cell_index = cell_index + nb[b];
						}
					}
				}// end of b
				// gamma_{b1}
				cell_index = 0;
				for (int b = 0; b < B; b++) {
					if (Drop_ind[b]) {
						//pro posal
						gamma_iter = -MCMC_Rng.rgamma(-10 * gamma[b][1], 0.1);
						//if(gamma_iter < 0){
//prior
						logr = (gamma_prior[0] - 1) * (log(gamma_iter / gamma[b][1])) + gamma_prior[1] * (gamma_iter - gamma[b][1]);
						//proposal
//numerator
						logr = logr - lgamma(-10 * gamma_iter) + (-10 * gamma_iter - 1) * log(-gamma[b][1]) - 10 * gamma_iter * log(10) + 10 * gamma[b][1];
						//denomerator
						logr = logr + lgamma(-10 * gamma[b][1]) - (-10 * gamma[b][1] - 1) * log(-gamma_iter) + 10 * gamma[b][1] * log(10) - 10 * gamma_iter;
						for (int i = 0; i < nb[b]; i++) {

#pragma omp parallel
							{
								// auxiliary variable to update gamma
								double temp_pres, temp_iter;
								double logr_thread = 0.0;
								int* X_thread, *Z_thread;
#pragma omp for
								for (int g = 0; g < G; g++) {
									X_thread = X[g];
									Z_thread = Z[g];
									//numerator
									temp_iter = gamma[b][0] + X_thread[cell_index] * gamma_iter;//prevent temp_iter is extremely large
									if (temp_iter > 0) {
										logr_thread += X_thread[cell_index] * gamma_iter * Z_thread[cell_index] - temp_iter - log(1 + exp(-temp_iter));
									}
									else {
										logr_thread += X_thread[cell_index] * gamma_iter * Z_thread[cell_index] - log(1 + exp(temp_iter));
									}
									//denomerator
									temp_pres = gamma[b][0] + X_thread[cell_index] * gamma[b][1];
									if (temp_pres > 0) {
										logr_thread += -X_thread[cell_index] * gamma[b][1] * Z_thread[cell_index] + temp_pres + log(1 + exp(-temp_pres));
									}
									else {
										logr_thread += -X_thread[cell_index] * gamma[b][1] * Z_thread[cell_index] + log(1 + exp(temp_pres));
									}
									//cout << " g = " << g << ", logr_thread = " << logr_thread << endl;
								}
#pragma omp critical
								{
									logr += logr_thread;

								}

							}
							cell_index++;
						}

						// cout << " logr = " << logr << endl;
						if (logr > log(MCMC_Rng.runif())) {
							gamma[b][1] = gamma_iter;
						}
					}
					else {
						cell_index = cell_index + nb[b];
					}

				}
			}
			////////////////////////////////////
			//  3) update alpha_g by MH  //
#pragma omp parallel for
			for (int g = 0; g < G; g++) {
				alpha[g] = _update_alpha(B, nb,//dimension
					mu_a[g], sigma_a,//prior
					W, alpha[g], beta[g], nu[g], delta, phi[g], //parameter
					X[g], MCMC_Rng);

			}
			////////////////////////
//  4) update L_gk  //
#pragma omp parallel for
			for (int g = 0; g < G; g++) {
				_update_l(K,//dimension
					p, tau0, tau1,//prior
					beta[g], MCMC_Rng,//parameter
					L[g]);
			}
			/////////////////////////////////////
//  5) update p and 6) update tau0  //
			if (IND_UPDATE_PTAU0 == 1) {
				// cout << "Update p and tau0." << endl;
				// auto start_pt = chrono::system_clock::now();
				int sum_L = 0;
				double sum_beta = 0.0;
#pragma omp parallel
				{
					int sum_L_thread = 0;
					double sum_beta_thread = 0.0;
					int *L_thread;
					double *beta_thread;
#pragma omp for
					for (int g = 0; g < G; g++) {
						L_thread = L[g];
						beta_thread = beta[g];
						for (int k = 1; k < K; k++) {
							sum_L_thread += L_thread[k];
							if (L_thread[k] == 0) {
								sum_beta_thread += pow(beta_thread[k], 2.0);
							}
						}
					}
#pragma omp critical
					{
						sum_L += sum_L_thread;
						sum_beta += sum_beta_thread;
					}
				}

				p_postdist[0] = p_prior[0] + sum_L;
				p_postdist[1] = p_prior[1] + G * (K - 1) - sum_L;
				p = MCMC_Rng.rbeta(p_postdist[0], p_postdist[1]);


				tau0_postdist[0] = tau0_prior[0] + (G * (K - 1) - sum_L) / 2.0;
				tau0_postdist[1] = tau0_prior[1] + sum_beta / 2.0;
				tau0 = 1.0 / MCMC_Rng.rgamma(tau0_postdist[0], 1.0 / tau0_postdist[1]);

			}
			///////////////////////
// 7) update beta  //
#pragma omp parallel for
			for (int g = 0; g < G; g++) {
				_update_beta(B, nb, K,//dimension 
					tau0, tau1, L[g],//prior
					W, alpha[g], nu[g], delta, phi[g], //parameter 	
					X[g], MCMC_Rng,//latent variable
					beta[g]);
			}
			///////////////////
			// 8) update nu  //
#pragma omp parallel for
			for (int g = 0; g < G; g++) {
				_update_nu(B, nb,
					mu_c, sigma_c,//prior
					W, alpha[g], beta[g], delta, phi[g],//parameter
					X[g], MCMC_Rng, //latent variable
					nu[g]);
			}
			/////////////////////
// 9) update delta /
			cell_index = 0;
			for (int b = 0; b < B; b++) {
				for (int i = 0; i < nb[b]; i++) {
					if (i > 0) {//let ind_n be the cell index
						delta_iter = MCMC_Rng.rnorm(delta[cell_index], 0.1);
						logr = 0.0;
						//prior
						logr += -pow(delta_iter - mu_d[cell_index], 2.0) / 2 / pow(sigma_d, 2.0);
						logr += pow(delta[cell_index] - mu_d[cell_index], 2.0) / 2 / pow(sigma_d, 2.0);
#pragma omp parallel
						{
							double logr_thread = 0.0;
							int *X_thread;
							double *beta_thread;
							double *nu_thread;
							double *phi_thread;
#pragma omp for
							for (int g = 0; g < G; g++) {
								X_thread = X[g];
								beta_thread = beta[g];
								nu_thread = nu[g];
								phi_thread = phi[g];
								int k = W[cell_index];
								//numerator
								logr_thread += delta_iter * X_thread[cell_index];
								logr_thread += -(phi_thread[b] + X_thread[cell_index]) * log(phi_thread[b] + exp(alpha[g] + beta_thread[k] + nu_thread[b] + delta_iter));
								//denomerator
								logr_thread += -delta[cell_index] * X_thread[cell_index];
								logr_thread += (phi_thread[b] + X_thread[cell_index]) * log(phi_thread[b] + exp(alpha[g] + beta_thread[k] + nu_thread[b] + delta[cell_index]));

							}
#pragma omp critical
							{
								logr += logr_thread;
							}
						}
						if (logr > log(MCMC_Rng.runif())) {
							delta[cell_index] = delta_iter;
						}
					}
					cell_index++;
				}
			}
#pragma omp parallel for
			for (int g = 0; g < G; g++) {
				_update_logmu(B, nb,
					W, alpha[g], beta[g], nu[g], delta,//parameter
					logmu[g]);
			}
			/////////////////////
// 10) update phi  //
#pragma omp parallel for
			for (int g = 0; g < G; g++) {
				_update_phi(B, nb,
					phi_prior,//prior
					logmu[g],//parameter
					X[g], MCMC_Rng, //latent variable
					phi[g]);
			}
			///////////////////
// 11) update w  //
			cell_index = 0;
			for (int b = 0; b < B; b++) {
				for (int i = 0; i < nb[b]; i++) {
					w_proposal = rand_cate(proposal_pi, MCMC_Rng);
					w_current = W[cell_index];

					if (w_proposal != w_current) {
						log_proposal = log(prop[b][w_proposal]);
						log_current = log(prop[b][w_current]);
						//calculate the posterior ratio in log scale
#pragma omp parallel
						{
							//double logr_thread = 0.0;
							double log_proposal_thread, log_current_thread;
							log_proposal_thread = 0.0;
							log_current_thread = 0.0;
							double temp_logmu;
							int *X_thread;
							double *beta_thread;
							double *nu_thread;
							double *phi_thread;
#pragma omp for
							for (int g = 0; g < G; g++) {
								X_thread = X[g];
								beta_thread = beta[g];
								nu_thread = nu[g];
								phi_thread = phi[g];

								temp_logmu = alpha[g] + beta_thread[w_proposal] + nu_thread[b] + delta[cell_index];
								log_proposal_thread += beta_thread[w_proposal] * X_thread[cell_index];
								log_proposal_thread += -(phi_thread[b] + X_thread[cell_index]) * log(phi_thread[b] + exp(temp_logmu));

								//ind_beta = j + w_current * _G;
								temp_logmu = alpha[g] + beta_thread[w_current] + nu_thread[b] + delta[cell_index];
								log_current_thread += beta_thread[w_current] * X_thread[cell_index];
								log_current_thread += -(phi_thread[b] + X_thread[cell_index]) * log(phi_thread[b] + exp(temp_logmu));

							}
#pragma omp critical
							{
								log_proposal += log_proposal_thread;
								log_current += log_current_thread;
							}
						}
						logr = log_proposal - log_current;
						if (logr > log(MCMC_Rng.runif())) {
							W[cell_index] = w_proposal;
						}
					}

					cell_index++;
				}
			}
#pragma omp parallel for
			for (int g = 0; g < G; g++) {
				_update_logmu(B, nb,
					W, alpha[g], beta[g], nu[g], delta,//parameter
					logmu[g]);
			}
			////////////////////////
//  12) update pi_bk  //
			cell_index = 0;
			for (int b = 0; b < B; b++) {
				for (int k = 0; k < K; k++) {
					count_w[k] = xi;
				}
				for (int i = 0; i < nb[b]; i++) {
					count_w[W[cell_index]] = count_w[W[cell_index]] + 1.0;
					cell_index++;
				}
				rand_Dir(count_w, K, MCMC_Rng, prop[b]);
			}
			if (iter == iter_noupdate - 1) {
				IND_UPDATE_PTAU0 = 1;
			}
			/////////////////////////////////////////
//  13) Record the posterior sampling  //
#pragma omp parallel for
			for (int g = 0; g < G; g++) {
				alpha_post[iter][g] = alpha[g];
			}
#pragma omp parallel
			{
				int q_thread;
#pragma omp for
				for (int g = 0; g < G; g++) {
					q_thread = g * K;
					for (int k = 0; k < K; k++) {
						beta_post[iter][q_thread] = beta[g][k];
						q_thread++;
					}
				}
			}
#pragma omp parallel
			{
				int q_thread;
#pragma omp for
				for (int g = 0; g < G; g++) {
					q_thread = g * B;
					for (int b = 0; b < B; b++) {
						nu_post[iter][q_thread] = nu[g][b];
						q_thread++;
					}
				}
			}
			for (int i = 0; i < N; i++) {
				delta_post[iter][i] = delta[i];
			}
			q = 0;
			for (int b = 0; b < B; b++) {
				gamma_post[iter][q] = gamma[b][0];
				q++;
				gamma_post[iter][q] = gamma[b][1];
				q++;
			}
#pragma omp parallel
			{
				int q_thread;
#pragma omp for
				for (int g = 0; g < G; g++) {
					q_thread = g * B;
					for (int b = 0; b < B; b++) {
						phi_post[iter][q_thread] = phi[g][b];
						q_thread++;
					}
				}
			}
			q = 0;
			for (int b = 0; b < B; b++) {
				for (int k = 0; k < K; k++) {
					pi_post[iter][q] = prop[b][k];
					q++;
				}
			}
			for (int i = 0; i < N; i++) {
				w_post[iter][i] = W[i];
			}
			p_post[iter] = p;
			tau0_post[iter] = tau0;
#pragma omp parallel
			{
				int q_thread;
#pragma omp for
				for (int g = 0; g < G; g++) {
					q_thread = g * K;
					for (int k = 0; k < K; k++) {
						l_post[iter][q_thread] = L[g][k];
						q_thread++;
					}
				}
			}
			// cout <<  "In " << iter << "-th iterations, the random number is " << MCMC_Rng.runif(0,1) << endl;

	              // output the progress bar
               if(iter % 100 == 0){
                       progress = tot_iter/iter_max;
                       int pos = barWidth * progress;
                       Rprintf("[");
                       for (int i = 0; i < barWidth; ++i) {
                           if (i < pos){
                                 Rprintf("=");
                           }else if(i == pos){
                                 Rprintf(">");
                           }else{
                                 Rprintf(" ");
                           }
               }
		       Rprintf("] Finish %.2fk/%.2fk iterations.\r", tot_iter/1000.0, iter_max/1000.0);
		       // fflush(stdout);
		       }
		       tot_iter = tot_iter + 1.0;

		}
		auto end_MCMC = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds_MCMC = end_MCMC - start_MCMC;
		// cout << "elapsed time of " << iter_num << " iterations of MCMC sampling for the " << t + 1 << "-th output is: " << elapsed_seconds_MCMC.count() << "s" << endl;
		///////////////////////////////////
		// output the posterior sampling //
		// cout << "Writing posterior sampling into the directory " << output_dir << endl;
		out_file = output_dir + "alpha_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			for (int g = 0; g < G; g++) {
				post_File << alpha_post[iter][g];
				post_File << " ";
			}
			post_File << endl;
		}
		post_File.close();

		out_file = output_dir + "beta_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			for (q = 0; q < G * K; q++) {
				post_File << beta_post[iter][q];
				post_File << " ";
			}
			post_File << endl;
		}
		post_File.close();

		out_file = output_dir + "nu_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			for (q = 0; q < G * B; q++) {
				post_File << nu_post[iter][q];
				post_File << " ";
			}
			post_File << endl;
		}
		post_File.close();
		out_file = output_dir + "delta_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			for (int i = 0; i < N; i++) {
				post_File << delta_post[iter][i];
				post_File << " ";

			}
			post_File << endl;
		}
		post_File.close();
		out_file = output_dir + "gamma_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			for (q = 0; q < B * 2; q++) {
				post_File << gamma_post[iter][q];
				post_File << " ";
			}
			post_File << endl;
		}
		post_File.close();
		out_file = output_dir + "phi_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			for (q = 0; q < G * B; q++) {
				post_File << phi_post[iter][q];
				post_File << " ";
			}
			post_File << endl;
		}
		post_File.close();
		out_file = output_dir + "pi_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			for (q = 0; q < B * K; q++) {
				post_File << pi_post[iter][q];
				post_File << " ";
			}
			post_File << endl;
		}
		post_File.close();
		out_file = output_dir + "w_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			for (int i = 0; i < N; i++) {
				post_File << w_post[iter][i];
				post_File << " ";
			}
			post_File << endl;
		}
		post_File.close();
		out_file = output_dir + "p_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			post_File << p_post[iter];
			post_File << endl;
		}
		post_File.close();
		out_file = output_dir + "tau0_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			post_File << tau0_post[iter];
			post_File << endl;
		}
		post_File.close();
		out_file = output_dir + "l_post.txt";
		post_File.open(out_file.c_str(), ios::out | ios::app);
		for (int iter = 0; iter < iter_num; iter++) {
			for (q = 0; q < G * K; q++) {
				post_File << l_post[iter][q];
				post_File << " ";
			}
			post_File << endl;
		}
		post_File.close();
		/*
		// output the imputed true read count
		if (t == out_times - 1) {
			out_file = output_dir + "x_imputed.txt";
			post_File.open(out_file.c_str(), ios::out | ios::app);
			for (int g = 0; g < G; g++) {
				for (int i = 0; i < N; i++) {
					post_File << X[g][i];
					post_File << " ";
				}
				post_File << endl;
			}
			post_File.close();
		}
		*/
	}
	
    progress = 1.0;
    int pos = barWidth * progress;
    Rprintf("[");
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos){
          Rprintf("=");
        }else if(i == pos){
          Rprintf(">");
        }else{
          Rprintf(" ");
        }
    }
    Rprintf("] Finish %.2fk/%.2fk iterations.\n\n", tot_iter/1000.0, iter_max/1000.0);
    // fflush(stdout);
	
	auto end_overall = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds_overall = end_overall - start_overall;
	// cout << "elapsed time of the overall algorithm is: " << elapsed_seconds_overall.count() << "s" << endl;
	//free the memory
	delete[] count_w;
	delete[] proposal_pi;
	delete[] p_postdist;
	delete[] tau0_postdist;
	for (int iter = 0; iter < iter_out; iter++) {
		delete[] alpha_post[iter];
	}
	delete[] alpha_post;
	for (int iter = 0; iter < iter_out; iter++) {
		delete[] beta_post[iter];
	}
	delete[] beta_post;
	for (int iter = 0; iter < iter_out; iter++) {
		delete[] nu_post[iter];
	}
	delete[] nu_post;
	for (int iter = 0; iter < iter_out; iter++) {
		delete[] delta_post[iter];
	}
	delete[] delta_post;
	for (int iter = 0; iter < iter_out; iter++) {
		delete[] gamma_post[iter];
	}
	delete[] gamma_post;
	for (int iter = 0; iter < iter_out; iter++) {
		delete[] phi_post[iter];
	}
	delete[] phi_post;
	for (int iter = 0; iter < iter_out; iter++) {
		delete[] pi_post[iter];
	}
	delete[] pi_post;
	for (int iter = 0; iter < iter_out; iter++) {
		delete w_post[iter];
	}
	delete[] w_post;
	delete[] p_post;
	delete[] tau0_post;
	for (int iter = 0; iter < iter_out; iter++) {
		delete[] l_post[iter];
	}
	delete[] l_post;
	for (int g = 0; g < G; g++) {
		delete[] Z[g];
	}
	delete[] Z;
	delete[] W;
	for (int g = 0; g < G; g++) {
		delete[] L[g];
	}
	delete[] L;
	for (int b = 0; b < B; b++) {
		delete[] prop[b];
	}
	delete[] prop;
	delete[] alpha;
	delete[] mu_a;
	for (int g = 0; g < G; g++) {
		delete[] beta[g];
	}
	delete[] beta;
	for (int g = 0; g < G; g++) {
		delete[] nu[g];
	}
	delete[] nu;
	delete[] mu_c;
	delete[] delta;
	delete[] mu_d;
	for (int g = 0; g < G; g++) {
		delete[] phi[g];
	}
	delete[] phi;
	for (int b = 0; b < B; b++) {
		delete[] gamma[b];
	}
	delete[] gamma;
}

void BUSseq_inference(int *Y_vec, int *Dim,
	int *nc, int *iter_infor,
	char **dir_output,
	double *fdr_threshold,
	// posterior mean, mode and standard deviation
	double *alpha_est, double *alpha_sd,
	double * beta_est, double *beta_sd,
	double *nu_est, double *nu_sd,
	double *delta_est, double *delta_sd,
	double *gamma_est, double * gamma_sd,
	double *phi_est, double *phi_sd,
	double *pi_est, double *pi_sd,
	double *tau0_est, double *tau0_sd,
	double *p_est, double *p_sd,
	int *w_est, double *PPI_est,
	int *D_est, double *BIC) {
	////////////////////////////////////////
	//  0. Deal with the input from R  //
	////////////////////////////////////////
	// Load the dimension information and raw count data
	int N = Dim[0];
	int G = Dim[1];
	int B = Dim[2];
	int K = Dim[3];
	int *nb = &(Dim[4]);
	int *di = &(Dim[4 + B]);

	bool *Drop_ind = new bool[B];
	bool All_Drop = true;
	for (int b = 0; b < B; b++) {
		if (di[b] == 1) {
			Drop_ind[b] = true;
		}
		else if (di[b] == 0) {
			Drop_ind[b] = false;
			All_Drop = false;
		}
		else {
		}
	}

	int **Y = new int*[G];
	for (int g = 0; g < G; g++) {
		Y[g] = &Y_vec[g * N];// for parallel G
	}

	// Set the number of cores for parallel
	omp_set_num_threads(nc[0]);

	// Set the number of iterations
	int iter_max = iter_infor[0]; // the overall iteration number 
	int iter_burnin = iter_infor[3];
	int n_iter = iter_max - iter_burnin;

	// Set the directory to output the posterior sampling
	string output_dir(dir_output[0]);
	output_dir = output_dir + "/";

#ifdef linux
int check = mkdir(output_dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif

#ifdef _WIN32
if (0 != access(output_dir.c_str(), 0))
    {
        // if this folder not exist, create a new one.
        mkdir(output_dir.c_str()); 
    }
#endif

	string out_file;

	// file variable to output the posterior sampling
	ofstream post_File;

	auto start_overall = chrono::system_clock::now();
	////////////////////////////
	// 1. Posterior inference //
	////////////////////////////
	////////////////////////////////////////////////
	//  1) load posterior sampling of parameters  //
	ifstream load_File;
	ofstream est_File;
	string load_name, est_name;
	double temp;
	int iter;
	string row_dropped;
	int q;

	// calculate alpha_est
	double *sum_alpha_sq = new double[G];
	for (int g = 0; g < G; g++) {
		sum_alpha_sq[g] = 0.0;
	}
	load_name = output_dir + "alpha_post.txt";
	load_File.open(load_name);
	iter = 0;

	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		for (int g = 0; g < G; g++) {
			load_File >> temp;
			alpha_est[g] += temp;
			sum_alpha_sq[g] += pow(temp, 2.0);
		}
	}
	for (int g = 0; g < G; g++) {
		alpha_est[g] = alpha_est[g] / n_iter;
		alpha_sd[g] = (sum_alpha_sq[g] - n_iter * pow(alpha_est[g], 2.0)) / (n_iter - 1);
	}
	load_File.close();

	// calculate beta_est
	double *sum_beta_sq = new double[K * G];

	for (q = 0; q < G * K; q++) {
		sum_beta_sq[q] = 0.0;
	}

	load_name = output_dir + "beta_post.txt";
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		for (q = 0; q < G * K; q++) {
			load_File >> temp;
			beta_est[q] += temp;
			sum_beta_sq[q] += pow(temp, 2.0);
		}
	}

	for (q = 0; q < G * K; q++) {
		beta_est[q] = beta_est[q] / n_iter;
		beta_sd[q] = (sum_beta_sq[q] - n_iter * pow(beta_est[q], 2.0)) / (n_iter - 1);
	}
	load_File.close();

	// calculate nu_est
	double *sum_nu_sq = new double[B * G];

	for (q = 0; q < G * B; q++) {
		sum_nu_sq[q] = 0.0;
	}

	load_name = output_dir + "nu_post.txt";
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		for (q = 0; q < G * B; q++) {
			load_File >> temp;
			nu_est[q] += temp;
			sum_nu_sq[q] += pow(temp, 2.0);
		}
	}

	for (q = 0; q < G * B; q++) {
		nu_est[q] = nu_est[q] / n_iter;
		nu_sd[q] = (sum_nu_sq[q] - n_iter * pow(nu_est[q], 2.0)) / (n_iter - 1);
	}
	load_File.close();

	// calculate delta_est
	double *sum_delta_sq = new double[N];
	for (int i = 0; i < N; i++) {
		sum_delta_sq[i] = 0.0;
	}
	load_name = output_dir + "delta_post.txt";
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		for (int i = 0; i < N; i++) {
			load_File >> temp;
			delta_est[i] += temp;
			sum_delta_sq[i] += pow(temp, 2.0);
		}
	}
	for (int i = 0; i < N; i++) {
		delta_est[i] = delta_est[i] / n_iter;
		delta_sd[i] = (sum_delta_sq[i] - n_iter * pow(delta_est[i], 2.0)) / (n_iter - 1);
	}
	load_File.close();

	// calculate gamma_est
	double *sum_gamma_sq = new double[B * 2];
	for (q = 0; q < 2 * B; q++) {
		sum_gamma_sq[q] = 0;
	}
	load_name = output_dir + "gamma_post.txt";
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		for (q = 0; q < B * 2; q++) {
			load_File >> temp;
			gamma_est[q] += temp;
			sum_gamma_sq[q] += pow(temp, 2.0);
		}
	}
	for (q = 0; q < B * 2; q++) {
		gamma_est[q] = gamma_est[q] / n_iter;
		gamma_sd[q] = (sum_gamma_sq[q] - n_iter * pow(gamma_est[q], 2.0)) / (n_iter - 1);
	}
	load_File.close();

	// calculate phi_est
	double *sum_phi_sq = new double[B * G];

	for (q = 0; q < G * B; q++) {
		sum_phi_sq[q] = 0.0;
	}

	load_name = output_dir + "phi_post.txt";
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		for (q = 0; q < G * B; q++) {
			load_File >> temp;
			phi_est[q] += temp;
			sum_phi_sq[q] += pow(temp, 2.0);
		}
	}

	for (q = 0; q < G * B; q++) {
		phi_est[q] = phi_est[q] / n_iter;
		phi_sd[q] = (sum_phi_sq[q] - n_iter * pow(phi_est[q], 2.0)) / (n_iter - 1);
	}
	load_File.close();

	// calculate pi_est
	double *sum_pi_sq = new double[B * K];

	for (q = 0; q < B * K; q++) {
		sum_pi_sq[q] = 0.0;
	}

	load_name = output_dir + "pi_post.txt";
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		for (q = 0; q < B * K; q++) {
			load_File >> temp;
			pi_est[q] += temp;
			sum_pi_sq[q] += pow(temp, 2.0);
		}
	}

	for (q = 0; q < B * K; q++) {
		pi_est[q] = pi_est[q] / n_iter;
		pi_sd[q] = (sum_pi_sq[q] - n_iter * pow(pi_est[q], 2.0)) / (n_iter - 1);
	}
	load_File.close();

	// calculate w_est
	int w_temp;

	int **count_w = new int*[N]; // count the frequency of cell types after burnin
	for (int i = 0; i < N; i++) {
		count_w[i] = new int[K];
		for (int k = 0; k < K; k++) {
			count_w[i][k] = 0;
		}
	}
	load_name = output_dir + "w_post.txt";
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		for (int i = 0; i < N; i++) {
			load_File >> w_temp;
			count_w[i][w_temp] ++;
		}
	}

	for (int i = 0; i < N; i++) {
		for (int k = 1; k < K; k++) {
			if (count_w[i][k] > count_w[i][w_est[i]]) {
				w_est[i] = k; // obtain the mode from burnin to the last iteration
			}
		}
	}
	load_File.close();

	// calculate p_est
	double sum_p_sq = 0.0;
	load_name = output_dir + "p_post.txt";
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		load_File >> temp;
		p_est[0] += temp;
		sum_p_sq += pow(temp, 2.0);

	}
	p_est[0] = p_est[0] / n_iter;
	p_sd[0] = (sum_p_sq - n_iter * pow(p_est[0], 2.0)) / (n_iter - 1);
	load_File.close();

	// calculate tau0_est
	double sum_tau0_sq = 0.0;
	load_name = output_dir + "tau0_post.txt";
	load_File.open(load_name);
	iter = 0;
	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		load_File >> temp;
		tau0_est[0] += temp;
		sum_tau0_sq += pow(temp, 2.0);

	}
	tau0_est[0] = tau0_est[0] / n_iter;
	tau0_sd[0] = (sum_tau0_sq - n_iter * pow(tau0_est[0], 2.0)) / (n_iter - 1);
		load_File.close();

	// calculate PPI_est
	load_name = output_dir + "l_post.txt";
	load_File.open(load_name);
	iter = 0;
	int temp_int;

	while (iter < iter_burnin) {
		getline(load_File, row_dropped);
		iter++;
	}
	for (; iter < iter_max; iter++) {
		for (q = 0; q < G * K; q++) {
			load_File >> temp_int;
			PPI_est[q] += temp_int;
		}
	}
	for (q = 0; q < G * K; q++) {
		PPI_est[q] = PPI_est[q] / n_iter;
	}
	load_File.close();

	// Extract intrinsic genes
	double postprob_DE_threshold;
	int count_gene = 0;
	postprob_DE_threshold = postprob_DE_thr_fun(PPI_est, fdr_threshold[0], G, K);

	for (int g = 0; g < G; g++) {
		D_est[g] = IG_index(PPI_est + g * K, postprob_DE_threshold, K);
		count_gene += D_est[g];
	}

	// Calculate likelihood and BIC

	omprng Loglike_Rng;
	Loglike_Rng.fixedSeed(12345);// for Monte carlo estimation of E(exp(gamma_0+gamma_1x)/(1+exp(gamma_0+gamma_1x)))
	int cell_index;
	double loglike_obs;
	// auxiliary variables
	double *lpy = new double[K]; // log(pi_bk prod_{g=1}^G Pr(Y_big = y_big | Theta ))
	double lpy_max, sum_lpy; // lr0_temp, sum_lr0, lpy_max, sum_lpy;
	// Set the number of cores for parallel
	omp_set_num_threads(nc[0]);

	// Calculate BIC 

	auto start_BIC = chrono::system_clock::now();
	loglike_obs = 0.0;
	cell_index = 0;

	for (int b = 0; b < B; b++) {
		if (Drop_ind[b]) {
			for (int i = 0; i < nb[b]; i++) {
				for (int k = 0; k < K; k++) {
					lpy[k] = log(pi_est[b * K + k]);
				}
				for (int k = 0; k < K; k++) {
					// auto start_bik = chrono::system_clock::now();
# pragma omp parallel
					{
						int read, x_max;
						double logmubikg, pbgk, logp, log1mp, lr0_temp, sum_lr0, lpy_thread;
						lpy_thread = 0.0;
# pragma omp for
						for (int g = 0; g < G; g++) {
							read = Y[g][cell_index];
							logmubikg = alpha_est[g] + beta_est[g * K + k] + nu_est[g * B + b] + delta_est[cell_index];
							pbgk = exp(logmubikg) / (exp(logmubikg) + phi_est[g * B+ b]);
							if (pbgk < exp(-100)) {
								logp = -100;
								log1mp = log(1 - pbgk);
							}
							else if (1 - pbgk < exp(-100)) {
								logp = log(pbgk);
								log1mp = -100;
							}
							else {
								logp = log(pbgk);
								log1mp = log(1 - pbgk);
							}

							if (read > 0) {
								lpy_thread += -log(1 + exp(gamma_est[b * 2] + gamma_est[b * 2 + 1] * read));
								lpy_thread += lgamma(phi_est[g * B + b] + read) - lgamma(read + 1) - lgamma(phi_est[g * B + b]);
								lpy_thread += read * logp + phi_est[g * B + b] * log1mp;
							}
							else {
								x_max = (int)3 * exp(logmubikg);
								lr0_temp = phi_est[g * B + b] * log1mp; //x=0
								sum_lr0 = lr0_temp;

								for (int x = 1; x < x_max; x++) {
									lr0_temp = gamma_est[b * 2] + gamma_est[b * 2 + 1] * x - log(1 + exp(gamma_est[b * 2] + gamma_est[b * 2 + 1] * x));
									lr0_temp += lgamma(phi_est[g * B + b] + x) - lgamma(x + 1) - lgamma(phi_est[g * B + b]);
									lr0_temp += x * logp + phi_est[g * B + b] * log1mp;
									if (lr0_temp > sum_lr0) {
										sum_lr0 = lr0_temp + log(1 + exp(sum_lr0 - lr0_temp));
									}
									else {
										sum_lr0 = sum_lr0 + log(1 + exp(lr0_temp - sum_lr0));
									}
								}
								lpy_thread += sum_lr0;
								//Rprintf("y %d %d %d = 0 and sum_lr0 is %f if the cell belongs to %d-th subtype.\n",b, i ,j ,sum_lr0, k);
							}
						}
# pragma omp critical
						{
							lpy[k] += lpy_thread;
						}
					}
				}// end of k 
				lpy_max = vec_max(lpy, K);
				sum_lpy = 0.0;
				for (int k = 0; k < K; k++) {
					sum_lpy = sum_lpy + exp(lpy[k] - lpy_max);
					//Rprintf("logproby[%d]=%f",k,lpy[k]);
				}
				loglike_obs += lpy_max + log(sum_lpy);
				cell_index++;
				// printf("Finish the %d-th cell, lpy_max= %f, sum_lpy = %f, loglike = %f.\n", cell_index, lpy_max, sum_lpy, loglike_obs);
			}// end of i
		}
		else {
			for (int i = 0; i < nb[b]; i++) {
				for (int k = 0; k < K; k++) {
					lpy[k] = log(pi_est[b * K + k]);
				}
				for (int k = 0; k < K; k++) {
					// auto start_bik = chrono::system_clock::now();
# pragma omp parallel
					{
						int read, x_max;
						double logmubikg, pbgk, logp, log1mp, lr0_temp, sum_lr0, lpy_thread;
						lpy_thread = 0.0;
# pragma omp for
						for (int g = 0; g < G; g++) {
							read = Y[g][cell_index];
							logmubikg = alpha_est[g] + beta_est[g * K + k] + nu_est[g * B + b] + delta_est[cell_index];
							pbgk = exp(logmubikg) / (exp(logmubikg) + phi_est[g * B + b]);
							if (pbgk < pow(0.1, 100)) {
								logp = -100;
								log1mp = log(1 - pbgk);
							}
							else if (1 - pbgk < pow(0.1, 100)) {
								logp = log(pbgk);
								log1mp = -100;
							}
							else {
								logp = log(pbgk);
								log1mp = log(1 - pbgk);
							}

							lpy_thread += lgamma(phi_est[g * B + b] + read) - lgamma(read + 1) - lgamma(phi_est[g * B+ b]);
							lpy_thread += read * logp + phi_est[g * B + b] * log1mp;
						}
# pragma omp critical
						{
							lpy[k] += lpy_thread;
						}
					}
				}// end of k 
				lpy_max = vec_max(lpy, K);
				sum_lpy = 0.0;
				for (int k = 0; k < K; k++) {
					sum_lpy = sum_lpy + exp(lpy[k] - lpy_max);
					//Rprintf("logproby[%d]=%f",k,lpy[k]);
				}
				loglike_obs += lpy_max + log(sum_lpy);
				cell_index++;
				// printf("Finish the %d-th cell, lpy_max= %f, sum_lpy = %f, loglike = %f.\n", cell_index, lpy_max, sum_lpy, loglike_obs);
			}// end of i
		} // end of else
	}// end of b

	int NumBatchDrop = 0;
	for (int b = 0; b < B; b++) {
		NumBatchDrop = NumBatchDrop + Drop_ind[b];
	}
	BIC[0] = -2.0 * loglike_obs + log(G * N) * ((B + G) * K + 2 * NumBatchDrop + G * (B * 2 - 1) + N - B);
	// all parameters of interest contain pi_{bk}, gamma_{b0(1)}, alpha_g, beta_{gk}, nu_{bg}, delta_{bi}, phi_{bg} 
	auto end_BIC = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds_BIC = end_BIC - start_BIC;

	delete[] sum_alpha_sq;
	delete[] sum_beta_sq;
	delete[] sum_nu_sq;
	delete[] sum_delta_sq;
	delete[] sum_gamma_sq;
	delete[] sum_phi_sq;
	delete[] sum_pi_sq;

	for (int i = 0; i < N; i++) {
		delete[] count_w[i];
	}
	delete[] count_w;
	delete[] lpy;
}

}


/////////////////////////
// Register routines //
/////////////////////////
// Describe the type and "style" of each argument
static R_NativePrimitiveArgType BUSseq_MCMC_t[] = {
    INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, STRSXP, REALSXP, INTSXP
};

static R_NativePrimitiveArgType BUSseq_inference_t[] = {
    INTSXP, INTSXP, INTSXP, INTSXP, STRSXP,  REALSXP, 
    REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, INTSXP, REALSXP,
    INTSXP, REALSXP
};

// Define R_CMethod
static const R_CMethodDef CEntries[] = {
  {"BUSseq_MCMC",          (DL_FUNC) &BUSseq_MCMC,   8,  BUSseq_MCMC_t},
  {"BUSseq_inference",          (DL_FUNC) &BUSseq_inference,   28,  BUSseq_inference_t},
  {NULL, NULL, 0}
};


// Register
void R_init_BUSseq(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries,NULL , NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
