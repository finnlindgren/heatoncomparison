#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <limits>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
using namespace std;

#include "kvpar.h"

#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <Rinternals.h>

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

extern "C" {
  extern void dgemm_(const char *transa, const char *transb,
		    const int *m, const int *n, const int *k,
		    const double *alpha, const double *a,
		    const int *lda, const double *b, const int *ldb,
		    const double *beta, double *c, const int *ldc);
   
  extern void  dcopy_(const int *n, const double *dx, const int *incx, double *dy, const int *incy);
  
  extern int dpotrf_(const char *uplo, int *n, double *a, int *lda, int *info);

  extern int dpotri_(const char *uplo, int *n, double *a, int *lda, int *info);

  extern void dsymm_(const char *side, const char *uplo, const int *m,
		     const int *n, const double *alpha,
		     const double *a, const int *lda,
		     const double *b, const int *ldb,
		     const double *beta, double *c, const int *ldc);

  extern double ddot_(int *n, double *x, int *incx, double *y, int *incy);

  extern void dgemv_(const char *trans, const int *m, const int *n, const double *alpha,
		     const double *a, const int *lda, const double *x, const int *incx,
		     const double *beta, double *y, const int *incy);
  
  extern void dsymv_(const char *uplo, const int *n, const double *alpha, const double *a, const int *lda,
		    const double *x, const int *incx, const double *beta, double *y, const int *incy);

  extern void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);

  extern void dtrmv_(const char *uplo, const char *transa, const char *diag, const int *n,
		     const double *a, const int *lda, double *b, const int *incx);

  extern void  dscal_(const int *n, double *da, double *dx, const int *incx);
  
}

void show(int *a, int n);
void show(double *a, int n);
void show(int *a, int r, int c);
void show(double *a, int r, int c);
void zeros(double *a, int n);
void writeRMatrix(string outfile, double * a, int nrow, int ncol);
void writeRMatrix(string outfile, int * a, int nrow, int ncol);
void mvrnorm(double *des, double *mu, double *cholCov, int dim);
double logit(double theta, double a, double b);
double logitInv(double z, double a, double b);
void dist(double *coords1, int n, double *coords2, int m, int p, double *d);
double rIG(double a, double b){
  return(1.0/rgamma(a, 1.0/b));
}

int main(int argc, char **argv){
  int h, i, j, k, l, s;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double negOne = -1.0;
  double zero = 0.0;
  char lower = 'L';
  char upper = 'U';
  char ntran = 'N';
  char ytran = 'T';
  char rside = 'R';
  char lside = 'L';
  
  string parfile;
  if(argc > 1)
    parfile = argv[1];
  else
    parfile = "pfile";
  
  kvpar par(parfile);
   
  //Get stuff
  int nThreads; par.getVal("n.threads", nThreads);
  int seed; par.getVal("seed", seed);
  int nSamples; par.getVal("n.samples", nSamples);
  int nReport; par.getVal("n.report", nReport);
  string outFile; par.getVal("out.file", outFile);
    
  vector<double> betaStarting; par.getVal("beta.starting", betaStarting);

  double tauSq_a = 2.0;
  double tauSq_b; par.getVal("tauSq.b", tauSq_b);
  double tauSqStarting; par.getVal("tauSq.starting", tauSqStarting);

  double sigmaSq_a = 2.0;
  double sigmaSq_b; par.getVal("sigmaSq.b", sigmaSq_b);
  double sigmaSqStarting; par.getVal("sigmaSq.starting", sigmaSqStarting);

  double phi_a; par.getVal("phi.a", phi_a);
  double phi_b; par.getVal("phi.b", phi_b);
  double phiStarting; par.getVal("phi.starting", phiStarting);
  double phiTuning; par.getVal("phi.tuning", phiTuning);
  
  omp_set_num_threads(nThreads);
  
  //set seed
  set_seed(123,seed);
  
  //n number of observations
  //p number of columns of X
  //m number of PP knots
  double *X = NULL; int n, p; par.getFile("X.file", X, n, p);
  double *y = NULL; par.getFile("y.file", y, n, i);
  double *coords = NULL; par.getFile("coords.file", coords, n, j);
  double *knots = NULL; int m; par.getFile("knot.file", knots, m, j);

  int nm = n*m;
  int mm = m*m;
  int pp = p*p;
  
  //parameters
  int nTheta = 3;
  int sigmaSqIndx = 0; 
  int tauSqIndx = 1; 
  int phiIndx = 2;

  //starting
  double *theta = new double[nTheta];
  theta[sigmaSqIndx] = sigmaSqStarting;
  theta[tauSqIndx] = tauSqStarting;
  theta[phiIndx] = phiStarting;
 
  //return stuff  
  double *betaSamples = new double[p*nSamples];
  for(i = 0; i < p; i++){
    betaSamples[i] = betaStarting[i];
  }
  double *wStrSamples = new double[m*nSamples];
  double *thetaSamples = new double[nTheta*nSamples];
 
  //dist and cov matrices
  double *D = new double[mm]; zeros(D, mm);
  double *C = new double[mm]; zeros(C, mm);
  double *d = new double[nm]; zeros(d, nm);
  double *c = new double[nm]; zeros(c, nm);
  
  dist(knots, m, knots, m, 2, D);
  dist(coords, n, knots, m, 2, d);
 
  //other stuff
  double logPostCand, logPostCurrent, logDetCurrent, logDetCand, accept = 0, status = 0;
  double *tmp_pp = new double[pp];
  double *tmp_p = new double[p];
  double *tmp_p2 = new double[p];
  double *tmp_m = new double[n];
  double *tmp_m2 = new double[n];
  double *tmp_n = new double[n];
  double *tmp_nm = new double[nm];
  double *tmp_mm = new double[mm];
  double *beta = new double[p];
  double *wStr = new double[m];
  double *w = new double[n];
  double tauSqInv, sigmaSqInv, phiCand;
  
  cout << "start sampling" << endl;

  //start Timers
  double *time = new double[2];
  double wall0 = get_wall_time();
  double cpu0 = get_cpu_time();

  double *xtx = new double[pp];
  dgemm_(&ytran, &ntran, &p, &p, &n, &one, X, &n, X, &n, &zero, xtx, &p);
 
  for(s = 0; s < nSamples; s++){
    
    ////////////////////////
    //update wStr and w
    ////////////////////////
    for(i = 0; i < mm; i++){
      C[i] = theta[sigmaSqIndx]*exp(-theta[phiIndx]*D[i]);
    }
    
    for(i = 0; i < nm; i++){
      c[i] = theta[sigmaSqIndx]*exp(-theta[phiIndx]*d[i]);
    }

    logDetCurrent = 0;
    dpotrf_(&lower, &m, C, &m, &info); if(info != 0){cout << "c++ error: dpotrf failed" << endl;}
    for(i = 0; i < m; i++){logDetCurrent += 2.0*log(C[i*m+i]);}
    dpotri_(&lower, &m, C, &m, &info); if(info != 0){cout << "c++ error: dpotri failed" << endl;}

    tauSqInv = 1.0/theta[tauSqIndx];

    dsymm_(&rside, &lower, &n, &m, &one, C, &m, c, &n, &zero, tmp_nm, &n);
    dgemm_(&ytran, &ntran, &m, &m, &n, &tauSqInv, tmp_nm, &n, tmp_nm, &n, &zero, tmp_mm, &m);

    daxpy_(&mm, &one, C, &inc, tmp_mm, &inc);

    dpotrf_(&lower, &m, tmp_mm, &m, &info); if(info != 0){cout << "c++ error: dpotrf failed" << endl;}
    dpotri_(&lower, &m, tmp_mm, &m, &info); if(info != 0){cout << "c++ error: dpotri failed" << endl;}

    dgemv_(&ntran, &n, &p, &one, X, &n, beta, &inc, &zero, tmp_n, &inc);
    for(i = 0; i < n; i++){
      tmp_n[i] =  y[i] - tmp_n[i];
    }

    dgemv_(&ytran, &n, &m, &tauSqInv, tmp_nm, &n, tmp_n, &inc, &zero, tmp_m, &inc);
   
    dsymv_(&lower, &m, &one, tmp_mm, &m, tmp_m, &inc, &zero, tmp_m2, &inc);
    
    dpotrf_(&lower, &m, tmp_mm, &m, &info); if(info != 0){cout << "c++ error: dpotrf failed" << endl;}
    mvrnorm(wStr, tmp_m2, tmp_mm, m);
    
    dsymv_(&lower, &m, &one, C, &m, wStr, &inc, &zero, tmp_m, &inc);
    dgemv_(&ntran, &n, &m, &one, c, &n, tmp_m, &inc, &zero, w, &inc);

    dcopy_(&m, wStr, &inc, &wStrSamples[s*m], &inc);
    
    ////////////////////////
    //update beta 
    ////////////////////////
    dcopy_(&pp, xtx, &inc, tmp_pp, &inc);
    dscal_(&pp, &tauSqInv, tmp_pp, &inc);
    
    dpotrf_(&lower, &p, tmp_pp, &p, &info); if(info != 0){cout << "c++ error: dpotrf failed" << endl;}
    dpotri_(&lower, &p, tmp_pp, &p, &info); if(info != 0){cout << "c++ error: dpotri failed" << endl;}

    for(i = 0; i < n; i++){
      tmp_n[i] = y[i] - w[i];
    }

    dgemv_(&ytran, &n, &p, &tauSqInv, X, &n, tmp_n, &inc, &zero, tmp_p, &inc);
    dsymv_(&lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc);
    
    dpotrf_(&lower, &p, tmp_pp, &p, &info); if(info != 0){cout << "c++ error: dpotrf failed" << endl;}
       
    mvrnorm(beta, tmp_p2, tmp_pp, p);

    dcopy_(&p, beta, &inc, &betaSamples[s*p], &inc);

    ////////////////////////
    //update vars
    ////////////////////////
    //sigmaSq
    dsymv_(&lower, &m, &theta[sigmaSqIndx], C, &m, wStr, &inc, &zero, tmp_m, &inc);
    
    theta[sigmaSqIndx] = rIG(sigmaSq_a + 0.5*m, sigmaSq_b + 0.5*ddot_(&m, wStr, &inc, tmp_m, &inc));
    
    //tauSq
    dgemv_(&ntran, &n, &p, &one, X, &n, beta, &inc, &zero, tmp_n, &inc);
    for(i = 0; i < n; i++){
      tmp_n[i] =  y[i] - tmp_n[i] - w[i];
    }

    theta[tauSqIndx] = rIG(tauSq_a + 0.5*n, tauSq_b + 0.5*ddot_(&n, tmp_n, &inc, tmp_n, &inc));
    
    ///////////////
    //update phi
    ///////////////
    //current
    dsymv_(&lower, &m, &one, C, &m, wStr, &inc, &zero, tmp_m, &inc); 
    logPostCurrent = -0.5*logDetCurrent - 0.5*ddot_(&m, wStr, &inc, tmp_m, &inc);
    logPostCurrent += log(theta[phiIndx] - phi_a) + log(phi_b - theta[phiIndx]); 

    //candidate
    phiCand = logitInv(rnorm(logit(theta[phiIndx], phi_a, phi_b), phiTuning), phi_a, phi_b);
    for(i = 0; i < mm; i++){
      C[i] = theta[sigmaSqIndx]*exp(-phiCand*D[i]);
    }
    
    logDetCand = 0;
    dpotrf_(&lower, &m, C, &m, &info); if(info != 0){cout << "c++ error: dpotrf failed" << endl;}
    for(i = 0; i < m; i++){logDetCand+= 2.0*log(C[i*m+i]);}
    dpotri_(&lower, &m, C, &m, &info); if(info != 0){cout << "c++ error: dpotri failed" << endl;}

    dsymv_(&lower, &m, &one, C, &m, wStr, &inc, &zero, tmp_m, &inc); 
    logPostCand = -0.5*logDetCand - 0.5*ddot_(&m, wStr, &inc, tmp_m, &inc);
    logPostCand += log(phiCand - phi_a) + log(phi_b - phiCand); 
    
    if(runif(0.0,1.0) <= exp(logPostCand - logPostCurrent)){
      theta[phiIndx] = phiCand;
      accept++;
    }
    
    dcopy_(&nTheta, theta, &inc, &thetaSamples[s*nTheta], &inc);

    if(status == nReport){
      cout << "percent complete: " << 100*s/nSamples << endl;
      cout << "metrop: " << 100.0*accept/s << endl;
      status = 0;
    }
    status++;

  }
  
  time[0] = get_wall_time() - wall0;
  time[1] = get_cpu_time() - cpu0;
  
  writeRMatrix(outFile+"-beta", betaSamples, p, nSamples);
  writeRMatrix(outFile+"-wStr", wStrSamples, m, nSamples);
  writeRMatrix(outFile+"-theta", thetaSamples, nTheta, nSamples);
  writeRMatrix(outFile+"-time", time, 2, 1);

  return(0);
}

void writeRMatrix(string outfile, double * a, int nrow, int ncol){

    ofstream file(outfile.c_str());
    if ( !file ) {
      cerr << "Data file could not be opened." << endl;
      exit(1);
    }
    
    for(int i = 0; i < nrow; i++){
      for(int j = 0; j < ncol-1; j++){
	file << setprecision(10) << fixed << a[j*nrow+i] << "\t";
      }
      file << setprecision(10) << fixed << a[(ncol-1)*nrow+i] << endl;    
    }
    file.close();
}


void writeRMatrix(string outfile, int* a, int nrow, int ncol){
  
  ofstream file(outfile.c_str());
  if ( !file ) {
    cerr << "Data file could not be opened." << endl;
    exit(1);
  }
  
  
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol-1; j++){
      file << fixed << a[j*nrow+i] << "\t";
    }
    file << fixed << a[(ncol-1)*nrow+i] << endl;    
  }
  file.close();
}

void mvrnorm(double *des, double *mu, double *cholCov, int dim){
  
  int i;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  
  for(i = 0; i < dim; i++){
    des[i] = rnorm(0, 1);
  }
 
  dtrmv_("L", "N", "N", &dim, cholCov, &dim, des, &inc);
  daxpy_(&dim, &one, mu, &inc, des, &inc);
}

void show(double *a, int n){
  for(int i = 0; i < n; i++)
    cout << setprecision(20) << fixed << a[i] << endl;
}


void show(int *a, int n){
  for(int i = 0; i < n; i++)
    cout << fixed << a[i] << endl;
}


void zeros(double *a, int n){
  for(int i = 0; i < n; i++)
    a[i] = 0.0;
}


void show(double *a, int r, int c){

  for(int i = 0; i < r; i++){
    for(int j = 0; j < c; j++){
      cout << fixed << a[j*r+i] << "\t";
    }
    cout << endl;
  }
}


void show(int *a, int r, int c){

  for(int i = 0; i < r; i++){
    for(int j = 0; j < c; j++){

      cout << fixed << a[j*r+i] << "\t";
    }
    cout << endl;
  }
}

double logit(double theta, double a, double b){
  return log((theta-a)/(b-theta));
}

double logitInv(double z, double a, double b){
  return b-(b-a)/(1+exp(z));
}

void dist(double *coords1, int n, double *coords2, int m, int p, double *d){

  double dist = 0.0;
  int i,j,k;

  for(i = 0; i < n; i++){
    for(j = 0; j < m; j++){
      for(k = 0; k < p; k++){
	d[j*n+i] += pow(coords1[n*k+i] - coords2[m*k+j], 2);
      }
      d[j*n+i] = sqrt(d[j*n+i]);
    }
  }
  
}

