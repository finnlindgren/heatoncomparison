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

#include "libs/kvpar.h"

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

}

static int rcmp(double x, double y, bool nalast)
{
    int nax = ISNAN(x), nay = ISNAN(y);
    if (nax && nay)	return 0;
    if (nax)		return nalast ? 1 : -1;
    if (nay)		return nalast ? -1 : 1;
    if (x < y)		return -1;
    if (x > y)		return 1;
    return 0;
}

void vsort_with_index(double *x, int *indx, int n)
{
    double v;
    int i, j, h, iv;

    for (h = 1; h <= n / 9; h = 3 * h + 1);
    for (; h > 0; h /= 3)
	for (i = h; i < n; i++) {
	    v = x[i]; iv = indx[i];
	    j = i;
	    while (j >= h && rcmp(x[j - h], v, TRUE) > 0)
		 { x[j] = x[j - h]; indx[j] = indx[j-h]; j -= h; }
	    x[j] = v; indx[j] = iv;
	}
}

void show(int *a, int n);
void show(double *a, int n);
void show(int *a, int r, int c);
void show(double *a, int r, int c);
void zeros(double *a, int n);
void writeRMatrix(string outfile, double * a, int nrow, int ncol);
void writeRMatrix(string outfile, int * a, int nrow, int ncol);
void mvrnorm(double *des, double *mu, double *cholCov, int dim);
void covTransInv(double *z, double *v, int m);
void covTrans(double *v, double *z, int m);
void covTrans(vector<double> v, double *z, int m);
void covTransInvExpand(double *v, double *z, int m);
void covExpand(double *v, double *z, int m);
double logit(double theta, double a, double b);
double logitInv(double z, double a, double b);
double dist2(double &a1, double &a2, double &b1, double &b2);

//Description: sorts vectors a and b of length n based on decreasing values of a.
void fSort(double *a, int *b, int n);

//Description: given a location's index i and number of neighbors m this function provides the index to i and number of neighbors in nnIndx
void getNNIndx(int i, int m, int &iNNIndx, int &iNN);

//Description: creates the nearest neighbor index given pre-ordered location coordinates.
//Input:
//n = number of locations
//m = number of nearest neighbors
//coords = ordered coordinates for the n locations
//Output:
//nnIndx = set of nearest neighbors for all n locations (on return)
//nnDist = euclidean distance corresponding to nnIndx (on return)
//nnIndxLU = nx2 look-up matrix with row values correspond to each location's index in nnIndx and number of neighbors (columns 1 and 2, respectively)
//Note: nnIndx and nnDist must be of length (1+m)/2*m+(n-m-1)*m on input. nnIndxLU must also be allocated on input.
void mkNNIndx(int n, int m, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU);

//void mkNNIndx2(int n, int m, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU);


//Description: using the fast mean-distance-ordered nn search by Ra and Kim 1993
//Input:
//ui = is the index for which we need the m nearest neighbors
//m = number of nearest neighbors
//n = number of observations, i.e., length of u
//sIndx = the NNGP ordering index of length n that is pre-sorted by u
//u = x+y vector of coordinates assumed sorted on input
//rSIndx = vector or pointer to a vector to store the resulting nn sIndx (this is at most length m for ui >= m)
//rNNDist = vector or point to a vector to store the resulting nn Euclidean distance (this is at most length m for ui >= m)  

double dmi(double *x, double *c, int inc){
    return pow(x[0]+x[inc]-c[0]-c[inc], 2);
}

double dei(double *x, double *c, int inc){
  return pow(x[0]-c[0],2)+pow(x[inc]-c[inc],2);
}

void fastNN(int m, int n, double *coords, int ui, double *u, int *sIndx, int *rSIndx, double *rSNNDist){
  
  int i,j,k;
  bool up, down;
  double dm, de;
  
  //rSNNDist will hold de (i.e., squared Euclidean distance) initially.
  for(i = 0; i < m; i++){
    rSNNDist[i] = std::numeric_limits<double>::infinity();
  }
  
  i = j = ui;
  
  up = down = true;
  
  while(up || down){
    
    if(i == 0){
      down = false;
    }

    if(j == (n-1)){
      up = false;
    }

    if(down){
      
      i--;
      
      dm = dmi(&coords[sIndx[ui]], &coords[sIndx[i]], n);
      
      if(dm > 2*rSNNDist[m-1]){
	down = false;
	
      }else{
	de = dei(&coords[sIndx[ui]], &coords[sIndx[i]], n);

	if(de < rSNNDist[m-1] && coords[sIndx[i]] < coords[sIndx[ui]]){
	  rSNNDist[m-1] = de;
	  rSIndx[m-1] = sIndx[i];
	  rsort_with_index(rSNNDist, rSIndx, m);
	}
	
      }
    }//end down
    
    if(up){
      
      j++;
      
      dm = dmi(&coords[sIndx[ui]], &coords[sIndx[j]], n);
      
      if(dm > 2*rSNNDist[m-1]){
	up = false;
	
      }else{
	de = dei(&coords[sIndx[ui]], &coords[sIndx[j]], n);

	if(de < rSNNDist[m-1] && coords[sIndx[j]] < coords[sIndx[ui]]){
	  rSNNDist[m-1] = de;
	  rSIndx[m-1] = sIndx[j];
	  rsort_with_index(rSNNDist, rSIndx, m);
	}
	
      }
      
    }//end up
    
  }
  
  for(i = 0; i < m; i++){
    rSNNDist[i] = sqrt(rSNNDist[i]);
  }


}


//Description: given the nnIndex this function fills uIndx for identifying those locations that have the i-th location as a neighbor.
//Input:
//n = number of locations
//m = number of nearest neighbors
//nnIndx = set of nearest neighbors for all n locations
//Output:
//uIndx = holds the indexes for locations that have each location as a neighbor
//uIndxLU = nx2 look-up matrix with row values correspond to each location's index in uIndx and number of neighbors (columns 1 and 2, respectively)
//Note: uIndx must be of length (1+m)/2*m+(n-m-1)*m on input. uINdxLU must also be allocated on input.
void mkUIndx(int n, int m, int* nnIndx, int* uIndx, int* uIndxLU);

//Description: writes nnIndex to file with each row corresponding to the ordered location coordinates (using R 1 offset). Each row corresponds the coordinate's index followed by its nearest neighbor indexes.
void writeRNNIndx(string outfile, int *nnIndx, int n, int m);

//Description: same as the other writeRNNIndx but this one uses the nnIndxLU look-up (duplicated function just for testing).
void writeRNNIndx(string outfile, int *nnIndx, int *nnIndxLU, int n);

//Description: same as writeRNNIndx but indexes in each row identify those locations that have the i-th row as a neighbor.
void writeRUIndx(string outfile, int *uIndx, int *uIndxLU, int n);

//Description: computes the quadratic term.
double Q(double *B, double *F, double *u, double *v, int n, int *nnIndx, int *nnIndxLU){
  
  double a, b, q = 0;
  int i, j;

#pragma omp parallel for private(a, b, j) reduction(+:q)
  for(i = 0; i < n; i++){
    a = 0;
    b = 0;
    for(j = 0; j < nnIndxLU[n+i]; j++){
      a += B[nnIndxLU[i]+j]*u[nnIndx[nnIndxLU[i]+j]];
      b += B[nnIndxLU[i]+j]*v[nnIndx[nnIndxLU[i]+j]];
    }
    q += (u[i] - a)*(v[i] - b)/F[i];
  }
  
  return(q);
}

//Description: update B and F.
double updateBF(double *B, double *F, double *c, double *C, double *D, double *d, int *nnIndxLU, int *CIndx, int n, double *theta, int tauSqIndx, int sigmaSqIndx, int phiIndx){
    
  int i, k, l;
  int info = 0;
  int inc = 1;
  double one = 1.0;
  double zero = 0.0;
  char lower = 'L';
  double logDet = 0;
    
#pragma omp parallel for private(k, l)
    for(i = 0; i < n; i++){
      if(i > 0){
	for(k = 0; k < nnIndxLU[n+i]; k++){
	  c[nnIndxLU[i]+k] = theta[sigmaSqIndx]*exp(-theta[phiIndx]*d[nnIndxLU[i]+k]);
	  for(l = 0; l <= k; l++){
	    C[CIndx[i]+l*nnIndxLU[n+i]+k] = theta[sigmaSqIndx]*exp(-theta[phiIndx]*D[CIndx[i]+l*nnIndxLU[n+i]+k]);
	    if(l == k){
	      C[CIndx[i]+l*nnIndxLU[n+i]+k] += theta[tauSqIndx];
	    }
	  }
	}
	dpotrf_(&lower, &nnIndxLU[n+i], &C[CIndx[i]], &nnIndxLU[n+i], &info); if(info != 0){cout << "c++ error: dpotrf failed" << endl;}
	dpotri_(&lower, &nnIndxLU[n+i], &C[CIndx[i]], &nnIndxLU[n+i], &info); if(info != 0){cout << "c++ error: dpotri failed" << endl;}
	dsymv_(&lower, &nnIndxLU[n+i], &one, &C[CIndx[i]], &nnIndxLU[n+i], &c[nnIndxLU[i]], &inc, &zero, &B[nnIndxLU[i]], &inc);
	F[i] = theta[sigmaSqIndx] - ddot_(&nnIndxLU[n+i], &B[nnIndxLU[i]], &inc, &c[nnIndxLU[i]], &inc) + theta[tauSqIndx];
      }else{
	B[i] = 0;
	F[i] = theta[sigmaSqIndx] + theta[tauSqIndx];
      }
    }
    
    for(i = 0; i < n; i++){
      logDet += 2*log(sqrt(F[i]));
    }

    return(logDet);
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
  
  bool debug = false;
  
  //Get stuff
  int nnIndexOnly; par.getVal("nn.index.only", nnIndexOnly);
  int nnFast; par.getVal("nn.fast", nnFast);
  int nThreads; par.getVal("n.threads", nThreads);
  int distThreads; par.getVal("dist.threads");
  int seed; par.getVal("seed", seed);
  int nSamples; par.getVal("n.samples", nSamples);
  int nReport; par.getVal("n.report", nReport);
  string outFile; par.getVal("out.file", outFile);
    
  vector<double> betaStarting; par.getVal("beta.starting", betaStarting);

  double tauSq_a = 2.0;
  double tauSq_b; par.getVal("tauSq.b", tauSq_b);
  double tauSqStarting; par.getVal("tauSq.starting", tauSqStarting);
  double tauSqTuning; par.getVal("tauSq.tuning", tauSqTuning);

  double sigmaSq_a = 2.0;
  double sigmaSq_b; par.getVal("sigmaSq.b", sigmaSq_b);
  double sigmaSqStarting; par.getVal("sigmaSq.starting", sigmaSqStarting);
  double sigmaSqTuning; par.getVal("sigmaSq.tuning", sigmaSqTuning);

  double phi_a; par.getVal("phi.a", phi_a);
  double phi_b; par.getVal("phi.b", phi_b);
  double phiStarting; par.getVal("phi.starting", phiStarting);
  double phiTuning; par.getVal("phi.tuning", phiTuning);

  omp_set_num_threads(distThreads);
  
  //set seed
  set_seed(123,seed);
  
  //m number of nearest neighbors
  //n number of observations
  //p number of columns of X
  
  int m; par.getVal("nn", m);
  
  //assuming coords, X, and y are already pre-ordered
  double *X = NULL; int n, p; par.getFile("X.file", X, n, p);
  double *y = NULL; par.getFile("y.file", y, n, i);
  double *coords = NULL; par.getFile("coords.file", coords, n, j);
  
  //allocated for the nearest neighbor index vector (note, first location has no neighbors).
  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  int *nnIndx = new int[nIndx];
  double *d = new double[nIndx];
  int *nnIndxLU = new int[2*n]; //first column holds the nnIndx index for the i-th location and the second columns holds the number of neighbors the i-th location has (the second column is a bit of a waste but will simplifying some parallelization).

  //start Timers
  double *time = new double[2];
  double wall0 = get_wall_time();
  double cpu0 = get_cpu_time();

  //make the neighbor index

  if(!nnFast){

    mkNNIndx(n, m, coords, nnIndx, d, nnIndxLU);
    //writeRNNIndx("nnIndx-slow", nnIndx, n, m);
    
  }else{
    
    int *sIndx = new int[n];
    double *u = new double[n];
    
    for(i = 0; i < n; i++){
      sIndx[i] = i;
      u[i] = coords[i]+coords[n+i];
    }
    
    cout << "sort" << endl;
    rsort_with_index(u, sIndx, n); 
    
    int iNNIndx, iNN;
    
    cout << "fastNN" << endl;
    //make nnIndxLU and fill nnIndx and d
#pragma omp parallel for private(iNNIndx, iNN)
    for(i = 0; i < n; i++){ //note this i indexes the u vector
      getNNIndx(sIndx[i], m, iNNIndx, iNN);
      nnIndxLU[sIndx[i]] = iNNIndx;
      nnIndxLU[n+sIndx[i]] = iNN;   
      fastNN(iNN, n, coords, i, u, sIndx, &nnIndx[iNNIndx], &d[iNNIndx]);
    } 
    //writeRNNIndx("nnIndx-fast", nnIndx, n, m);
  }
  
  if(nnIndexOnly){
    exit(1);
  }

  //parameters
  int nTheta = 3;
  int sigmaSqIndx = 0; 
  int tauSqIndx = 1; 
  int phiIndx = 2;

  //starting
  double *theta = new double[nTheta];
  double *thetaCand = new double[nTheta];
  theta[sigmaSqIndx] = sigmaSqStarting;
  theta[tauSqIndx] = tauSqStarting;
  theta[phiIndx] = phiStarting;
 
  //tuning
  double *tuning = new double[nTheta];
  tuning[sigmaSqIndx] = sigmaSqTuning;
  tuning[tauSqIndx] = tauSqTuning;
  tuning[phiIndx] = phiTuning;
  
  //return stuff  
  double *betaSamples = new double[p*nSamples];
  double *thetaSamples = new double[nTheta*nSamples];
 
  //other stuff
  double *B = new double[nIndx];
  double *F = new double[n];
  double *c = new double[nIndx];

  int *CIndx = new int[2*n]; //index for D and C.
  for(i = 0, j = 0; i < n; i++){//zero should never be accessed
    j += nnIndxLU[n+i]*nnIndxLU[n+i];
    if(i == 0){
      CIndx[n+i] = 0;
      CIndx[i] = 0;
    }else{
      CIndx[n+i] = nnIndxLU[n+i]*nnIndxLU[n+i]; 
      CIndx[i] = CIndx[n+i-1] + CIndx[i-1];
    }
  }
 
  double *C = new double[j]; zeros(C, j);
  double *D = new double[j]; zeros(D, j);

  for(i = 0; i < n; i++){
    for(k = 0; k < nnIndxLU[n+i]; k++){   
      for(l = 0; l <= k; l++){
  	D[CIndx[i]+l*nnIndxLU[n+i]+k] = dist2(coords[nnIndx[nnIndxLU[i]+k]], coords[n+nnIndx[nnIndxLU[i]+k]], coords[nnIndx[nnIndxLU[i]+l]], coords[n+nnIndx[nnIndxLU[i]+l]]);
      }
    }
  }
  
  if(debug){
    for(i = 0; i < n; i++){
      show(&d[nnIndxLU[i]], 1, nnIndxLU[n+i]);
      cout << "--" << endl;
      show(&D[CIndx[i]], nnIndxLU[n+i], nnIndxLU[n+i]);
      cout << "-------" << endl;
    }
  }

  //other stuff
  double logPostCand, logPostCurrent, logDetCurrent, logDetCand, QCurrent, QCand, accept = 0, status = 0;
  int pp = p*p;
  double *tmp_pp = new double[pp];
  double *tmp_p = new double[p];
  double *tmp_p2 = new double[p];
  double *beta = new double[p];
  double *tmp_n = new double[n];

  bool thetaUpdate = true;

  //update B and F
  logDetCurrent = updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, theta, tauSqIndx, sigmaSqIndx, phiIndx);

  for(i = 0; i < p; i++){betaSamples[i] = betaStarting[i];}
  dgemv_(&ntran, &n, &p, &one, X, &n, betaSamples, &inc, &zero, tmp_n, &inc);
  daxpy_(&n, &negOne, y, &inc, tmp_n, &inc);
  QCurrent = Q(B, F, tmp_n, tmp_n, n, nnIndx, nnIndxLU);

  omp_set_num_threads(nThreads);

  cout << "start sampling" << endl;

  for(s = 0; s < nSamples; s++){
    
    if(thetaUpdate){
      
      thetaUpdate = false;
      
      ///////////////
      //update beta 
      ///////////////
      for(i = 0; i < p; i++){
	tmp_p[i] = Q(B, F, &X[n*i], y, n, nnIndx, nnIndxLU);
	for(j = 0; j <= i; j++){
	  tmp_pp[j*p+i] = Q(B, F, &X[n*j], &X[n*i], n, nnIndx, nnIndxLU);
	}
      }
      
      dpotrf_(&lower, &p, tmp_pp, &p, &info); if(info != 0){cout << "c++ error: dpotrf failed" << endl;}
      dpotri_(&lower, &p, tmp_pp, &p, &info); if(info != 0){cout << "c++ error: dpotri failed" << endl;}
      
      dsymv_(&lower, &p, &one, tmp_pp, &p, tmp_p, &inc, &zero, tmp_p2, &inc);
      
      dpotrf_(&lower, &p, tmp_pp, &p, &info); if(info != 0){cout << "c++ error: dpotrf failed" << endl;}
      
    }
    
    mvrnorm(&betaSamples[s*p], tmp_p2, tmp_pp, p);
      
    ///////////////
    //update theta
    ///////////////
    dgemv_(&ntran, &n, &p, &one, X, &n, &betaSamples[s*p], &inc, &zero, tmp_n, &inc);
    daxpy_(&n, &negOne, y, &inc, tmp_n, &inc);
    
    //current    
    logDetCurrent = updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, theta, tauSqIndx, sigmaSqIndx, phiIndx);

    QCurrent = Q(B, F, tmp_n, tmp_n, n, nnIndx, nnIndxLU);

    logPostCurrent = -0.5*logDetCurrent - 0.5*QCurrent;
    logPostCurrent += log(theta[phiIndx] - phi_a) + log(phi_b - theta[phiIndx]); 
    logPostCurrent += -1.0*(1.0+sigmaSq_a)*log(theta[sigmaSqIndx])-sigmaSq_b/theta[sigmaSqIndx]+log(theta[sigmaSqIndx]);
    logPostCurrent += -1.0*(1.0+tauSq_a)*log(theta[tauSqIndx])-tauSq_b/theta[tauSqIndx]+log(theta[tauSqIndx]);
    
    //candidate
    thetaCand[phiIndx] = logitInv(rnorm(logit(theta[phiIndx], phi_a, phi_b), phiTuning), phi_a, phi_b);
    thetaCand[sigmaSqIndx] = exp(rnorm(log(theta[sigmaSqIndx]), sigmaSqTuning));
    thetaCand[tauSqIndx] = exp(rnorm(log(theta[tauSqIndx]), tauSqTuning));

    //update B and F
    logDetCand = updateBF(B, F, c, C, D, d, nnIndxLU, CIndx, n, thetaCand, tauSqIndx, sigmaSqIndx, phiIndx);

    QCand = Q(B, F, tmp_n, tmp_n, n, nnIndx, nnIndxLU);

    logPostCand = -0.5*logDetCand - 0.5*QCand;
    logPostCand += log(thetaCand[phiIndx] - phi_a) + log(phi_b - thetaCand[phiIndx]); 
    logPostCand += -1.0*(1.0+sigmaSq_a)*log(thetaCand[sigmaSqIndx])-sigmaSq_b/thetaCand[sigmaSqIndx]+log(thetaCand[sigmaSqIndx]);
    logPostCand += -1.0*(1.0+tauSq_a)*log(thetaCand[tauSqIndx])-tauSq_b/thetaCand[tauSqIndx]+log(thetaCand[tauSqIndx]);
    
    if(runif(0.0,1.0) <= exp(logPostCand - logPostCurrent)){
      thetaUpdate = true;
      dcopy_(&nTheta, thetaCand, &inc, theta, &inc);
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

void covTransInv(double *z, double *v, int m){
  int i, j, k;

  for(i = 0, k = 0; i < m; i++){
    for(j = i; j < m; j++, k++){
      v[k] = z[k];
      if(i == j)
	v[k] = exp(z[k]);
    }
  }

}

void covTrans(double *v, double *z, int m){
  int i, j, k;

  for(i = 0, k = 0; i < m; i++){
    for(j = i; j < m; j++, k++){
      z[k] = v[k];
      if(i == j)
	z[k] = log(v[k]);
    }
  }

}

void covTrans(vector<double> v, double *z, int m){
  int i, j, k;

  for(i = 0, k = 0; i < m; i++){
    for(j = i; j < m; j++, k++){
      z[k] = v[k];
      if(i == j)
	z[k] = log(v[k]);
    }
  }

}

void covTransInvExpand(double *v, double *z, int m){
  int i, j, k;
  
  zeros(z, m*m);
  for(i = 0, k = 0; i < m; i++){
    for(j = i; j < m; j++, k++){
      z[i*m+j] = v[k];
      if(i == j)
	z[i*m+j] = exp(z[i*m+j]);
    }
  }
  
}

void covExpand(double *v, double *z, int m){
  int i, j, k;
  
  zeros(z, m*m);
  for(i = 0, k = 0; i < m; i++){
    for(j = i; j < m; j++, k++){
      z[i*m+j] = v[k];
    }
  }
  
}

double logit(double theta, double a, double b){
  return log((theta-a)/(b-theta));
}

double logitInv(double z, double a, double b){
  return b-(b-a)/(1+exp(z));
}

double dist2(double &a1, double &a2, double &b1, double &b2){
  return(sqrt(pow(a1-b1,2)+pow(a2-b2,2)));
}

void fSort(double *a, int *b, int n){
  
  int j, k, l;
  double v;
  
  for(j = 1; j <= n-1; j++){
    k = j;  
    while(k > 0 && a[k] < a[k-1]) {
      v = a[k]; l = b[k];
      a[k] = a[k-1]; b[k] = b[k-1];
      a[k-1] = v; b[k-1] = l;
      k--;
    }
  }
}



void getNNIndx(int i, int m, int &iNNIndx, int &iNN){
  
  if(i == 0){
    iNNIndx = 0;//this should never be accessed
    iNN = 0;
    return;
  }else if(i < m){
    iNNIndx = static_cast<int>(static_cast<double>(i)/2*(i-1));
    iNN = i;
    return;
  }else{
    iNNIndx = static_cast<int>(static_cast<double>(m)/2*(m-1)+(i-m)*m);
    iNN = m;
    return;
  } 
}


void mkNNIndx(int n, int m, double *coords, int *nnIndx, double *nnDist, int *nnIndxLU){
  
  int i, j, iNNIndx, iNN;
  double d;
  
  int nIndx = static_cast<int>(static_cast<double>(1+m)/2*m+(n-m-1)*m);
  
  for(i = 0; i < nIndx; i++){
    nnDist[i] = std::numeric_limits<double>::infinity();
  }
  
  #pragma omp parallel for private(j, iNNIndx, iNN, d)
  for(i = 0; i < n; i++){ 
    getNNIndx(i, m, iNNIndx, iNN);
    nnIndxLU[i] = iNNIndx;
    nnIndxLU[n+i] = iNN;   
    if(i != 0){  
      for(j = 0; j < i; j++){	
	d = dist2(coords[i], coords[n+i], coords[j], coords[n+j]);	
	if(d < nnDist[iNNIndx+iNN-1]){	  
	  nnDist[iNNIndx+iNN-1] = d;
	  nnIndx[iNNIndx+iNN-1] = j;
	  fSort(&nnDist[iNNIndx], &nnIndx[iNNIndx], iNN); 	  
	}	
      }
    } 
  }
  
}


void writeRNNIndx(string outfile, int *nnIndx, int n, int m){
  
  ofstream file(outfile.c_str());
  if(!file){
    cerr << "Data file could not be opened." << endl;
    exit(1);
  }
  
  int i, j, a, b;
  
  for(i = 0; i < n; i++){
    if(i != 0){    
      getNNIndx(i, m, a, b);    
      file << i+1 << " ";
      for(j = 0; j < b; j++){
	if(j+1 == b){
	  file << nnIndx[a+j]+1;
	}else{
	  file << nnIndx[a+j]+1 << ",";
	}
      }
      file << endl;
    }
  }
  file.close();
}

void writeRNNIndx(string outfile, int *nnIndx, int *nnIndxLU, int n){
  
  ofstream file(outfile.c_str());
  if(!file){
    cerr << "Data file could not be opened." << endl;
    exit(1);
  }
 
  int i, j;
  
  for(i = 0; i < n; i++){
    if(nnIndxLU[n+i] > 0){//i.e., not i = 0
      file << i+1 << " ";
      for(j = 0; j < nnIndxLU[n+i]; j++){
	if(j+1 == nnIndxLU[n+i]){
	  file << nnIndx[nnIndxLU[i]+j]+1;
	}else{	
	  file << nnIndx[nnIndxLU[i]+j]+1 << ",";	
	}
      }
      file << endl;
    }
  }
  
  file.close();
}

void writeRUIndx(string outfile, int *uIndx, int *uIndxLU, int n){
  
  ofstream file(outfile.c_str());
  if(!file){
    cerr << "Data file could not be opened." << endl;
    exit(1);
  }

  int i, j;
  
  for(i = 0; i < n; i++){
    file << i+1 << " ";
    for(j = 0; j < uIndxLU[n+i]; j++){
      if(j+1 == uIndxLU[n+i]){
	file << uIndx[uIndxLU[i]+j]+1;
      }else{	
	file << uIndx[uIndxLU[i]+j]+1 << ",";	
      }
    }
    file << endl;
  }
  
  file.close();
}

void mkUIndx(int n, int m, int* nnIndx, int* uIndx, int* uIndxLU){ 
  
  int iNNIndx, iNN, i, j, k, l, h;
  
  for(i = 0, l = 0; i < n; i++){    
    uIndxLU[i] = l; 
    for(j = 0, h = 0; j < n; j++){   
      getNNIndx(j, m, iNNIndx, iNN);  
      for(k = 0; k < iNN; k++){      	
	if(nnIndx[iNNIndx+k] == i){
	  uIndx[l+h] = j;
	  h++;
	}    
      }
    }
    l += h;
    uIndxLU[n+i] = h; 
  }
}

