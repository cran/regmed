#include <RcppArmadillo.h>

// for print_matrix
#include <iostream>
#include <iomanip>

// #include <cmath>

using namespace std;
using namespace Rcpp;
using namespace arma;

///////////////  Function Prototypes //////////////////////////////////////////
double soft_threshold(double z, double lambda); 
// note that soft_threshold function is defined in regmed.cpp

bool isNaN(double x) { 
  return x != x;
}

arma::mat blockDiag(arma::mat &A, arma::mat &B, arma::mat &C);

arma::mat col_outerprod_row(arma::mat &A, arma::mat &B, int a_col, int b_row);

arma::mat compute_E(arma::mat &Alpha, arma::mat &Beta, arma::mat &Delta,
        arma::mat &Varx, arma::mat &Varm, arma::mat &Vary);

arma::mat compute_B(arma::mat &Alpha, arma::mat &Beta, arma::mat &Delta);

arma::mat compute_ImA(arma::mat &Alpha, arma::mat &Beta, arma::mat &Delta);

double logdet_Var(arma::mat Var);

arma::mat inverse_ImpCov(arma::mat &Alpha, arma::mat &Beta, arma::mat &Delta,
        arma::mat &Varx_inv, arma::mat &Varm_inv, arma::mat &Vary_inv);

double grad_Amat(int index_i,
        int index_j,
        arma::mat & Alpha,
        arma::mat & Beta,
        arma::mat & Delta,
        arma::mat & SampCov,
        arma::mat & S,
        arma::mat & Varx_inv, 
        arma::mat & Varm_inv,
        arma::mat & Vary_inv);

double grad_AmatV2(int index_i,
        int index_j,
        arma::mat & Alpha,
        arma::mat & Beta,
        arma::mat & Delta,
        arma::mat & SampCov,
        arma::mat & S,
        arma::mat & Varx_inv, 
        arma::mat & Varm_inv,
        arma::mat & Vary_inv,
        arma::mat & Varx,
        arma::mat & Varm,
        arma::mat & Vary);

double grad_Smat(int index_i,
        int index_j,
        arma::mat & Alpha,
        arma::mat & Beta,
        arma::mat & Delta,
        arma::mat & SampCov,
        arma::mat & Varx_inv,
        arma::mat & Varm_inv,
        arma::mat & Vary_inv);

double penalty(arma::mat &Alpha,
        arma::mat &Beta,
        arma::mat &Delta,
        double lambda);

double update_param(double gradient, double param, double step,
        double lambda);

int count_df(arma::mat &A);

int count_df_Vary(arma::mat &V);

arma::mat gchol(arma::mat matrix);

arma::mat gchol_inv(arma::mat matrix);

void print_vec(arma::vec x);

void print_mat(arma::mat &mat);


///////////////////////  Functions ///////////////////////////////////////////

arma::mat blockDiag(arma::mat &A, arma::mat &B, arma::mat &C) {

    // make block diag matrix as
    //   A  0  0
    //   0  B  0
    //   0  0  C

    // A, B, C are square matrices

    int nA = A.n_rows;
    int nB = B.n_rows;
    int nC = C.n_rows;
    int nr = nA + nB + nC;

    arma::mat D(nr, nr, fill::zeros);

    // X.submat( first_row, first_col, last_row, last_col )

    D.submat(0, 0, nA - 1, nA - 1) = A;
    D.submat(nA, nA, nA + nB - 1, nA + nB - 1) = B;
    D.submat(nA + nB, nA + nB, nA + nB + nC - 1, nA + nB + nC - 1) = C;

    return D;
}

arma::mat col_outerprod_row(arma::mat & A, arma::mat & B, int a_col, int b_row) {

    // For matrices A and B, select column a_col of A (call this a) and
    // select row b_row of B (call this b) and compute outer product matrix
    // C = a (x) b

    int nrowA = A.n_rows;
    int ncolA = A.n_cols;

    int nrowB = B.n_rows;
    int ncolB = B.n_cols;

    arma::mat C(nrowA, ncolB, fill::zeros);

    for (int r = 0; r < nrowA; r++) {
        for (int c = 0; c < ncolB; c++) {
            C(r, c) = A(r, a_col) * B(b_row, c);
        }
    }
    return C;
}

arma::mat compute_E(arma::mat &Alpha, arma::mat &Beta, arma::mat &Delta,
        arma::mat &Varx, arma::mat &Varm, arma::mat &Vary) {

    // Compared with efficiency for computing E = B * S * B.t(), the 
    // implied covar matrix, and this version is about 3 times faster

    // Alpha is nm x nx for mediators by x's
    // Beta  is nt x nm for traits by mediators
    // Delta is nt x nx for traits by x's
    
    // E is symmetric matrix, with lower blocks shown below
    // E rows: x, mediators, traits
    //         0:nx-1, nx:nx+nm-1, nx+nm:nx+nm+nt-1

    // E =  Varx       |
    //      Alpha*Varx | Alpha*Varx*Alpha' + Varm      |
    //      Gamma*Varx | Gamma*Varx*Alpha' + Beta*Varm | Gamma*Varx*Gamma' + Beta*Varm*Beta' + Vary

    arma::mat Gamma = Delta + Beta * Alpha; // dim = nt x nx
    arma::mat AlphaVarx = Alpha * Varx;     // dim = nm x nx
    arma::mat GammaVarx = Gamma * Varx;     // dim = nt x nx
    arma::mat BetaVarm  = Beta * Varm;       // dim = nt x nm

    int nx = Varx.n_rows;
    int nm = Varm.n_rows;
    int nt = Vary.n_rows;
    
    int nrowE = nx + nm + nt;
    arma::mat E(nrowE, nrowE, fill::zeros);

    // syntax: X.submat( first_row, first_col, last_row, last_col )

    E.submat(0, 0, nx - 1, nx - 1) = Varx;
    E.submat(nx, 0, nx + nm - 1, nx - 1) = AlphaVarx;
    E.submat(nx + nm, 0, nx + nm + nt - 1, nx - 1) = GammaVarx;

    E.submat(nx, nx, nx + nm - 1, nx + nm - 1) = AlphaVarx * Alpha.t() + Varm;
    E.submat(nx + nm, nx, nx + nm + nt - 1, nx + nm - 1) = GammaVarx * Alpha.t() + BetaVarm;

    E.submat(nx + nm, nx + nm, nx + nm + nt - 1, nx + nm + nt - 1) = GammaVarx * Gamma.t() + BetaVarm * Beta.t() + Vary; 
    
    // Fill in upper triangle
    for (int i = 0; i < (nrowE - 1); i++) {
        for (int j = (i + 1); j < nrowE; j++) {
            E(i, j) = E(j, i);
        }
    }

    return E;

}

arma::mat compute_B(arma::mat &Alpha, arma::mat &Beta, arma::mat &Delta) {

    // compute  B = inv(I-A), but compute B with analytic solution
    // Alpha is nm x nx (nm = no. mediators; nx = no. snps)
    // Delta is  nt x nx (nt = no. traits)
    // Beta is nt x nm

    // B =     I |  0    | 0
    //     Alpha |  I    | 0
    //     Gamma | Beta  | I

    // where Gamma = Delta + Beta * Alpha

    arma::mat Gamma = Delta + Beta * Alpha;

    int nm = Alpha.n_rows;
    int nx = Alpha.n_cols;
    int nt = Delta.n_rows;
    int ntot = nx + nm + nt;

    arma::mat B(ntot, ntot);
    B.eye();

    // syntax: X.submat( first_row, first_col, last_row, last_col )

    B.submat(nx, 0, nx + nm - 1, nx - 1) = Alpha;

    B.submat(nx + nm, 0, nx + nm + nt - 1, nx - 1) = Gamma;

    B.submat(nx + nm, nx, nx + nm + nt - 1, nx + nm - 1) = Beta;


    return B;
}

arma::mat compute_ImA(arma::mat &Alpha, arma::mat &Beta, arma::mat &Delta) {

    // compute I-A
 
    // I-A=     I |  0    | 0
    //     -Alpha |  I    | 0
    //     -Delta | -Beta | I


    int nm = Alpha.n_rows;
    int nx = Alpha.n_cols;
    int nt = Delta.n_rows;
    int ntot = nx + nm + nt;

    arma::mat ImA(ntot, ntot);
    ImA.eye();

    // syntax: X.submat( first_row, first_col, last_row, last_col )

    ImA.submat(nx, 0, nx + nm - 1, nx - 1) = -Alpha;

    ImA.submat(nx + nm, 0, nx + nm + nt - 1, nx - 1) = -Delta;

    ImA.submat(nx + nm, nx, nx + nm + nt - 1, nx + nm - 1) = -Beta;

    return ImA;
}

double logdet_Var(arma::mat Var) {

    // compute logdet of square symmetric matrix

    // old with Armadillo chol: arma::mat lower = chol(Var, "lower");
    // new with gchol
    arma::mat lower = gchol(Var);
    double eps = 0.00001;
    int nr = Var.n_rows;
    double logdet = 0.0;
    for (int i = 0; i < nr; i++) {
        if (lower(i, i) > eps) {
            logdet += log(lower(i, i));
        }
    }
    logdet = 2.0 * logdet;
    return logdet;
}

arma::mat inverse_ImpCov(arma::mat &Alpha, arma::mat &Beta, arma::mat &Delta,
        arma::mat &Varx_inv, arma::mat &Varm_inv, arma::mat &Vary_inv) {
    
    // This version is about 5-6 times faster than simpler version using
    // ImpCov_inv = ImA.t() * Sinv * ImA;


    // inverse of implied covar matrix = (I-A)' Sinv (I-A)
    //
    // (I-A) matrix
    //                I      0   0
    //           -Alpha      I   0
    //           -Delta  -Beta   I

    // Sinv matrix is block diagonal
    //   Varx_inv
    //             Varm_inv
    //                        Vary_inv

    int nx = Alpha.n_cols;
    int nm = Alpha.n_rows;
    int nt = Delta.n_rows;
    int ntot = nx + nm + nt;

    arma::mat ImpCov_inv(ntot, ntot, fill::zeros);

    arma::mat VmAlpha = Varm_inv * Alpha;
    arma::mat VyBeta  = Vary_inv * Beta;
    arma::mat VyDelta = Vary_inv * Delta;

    // X.submat( first_row, first_col, last_row, last_col )

    // Row and Col blocks
    // Block 1,1
    ImpCov_inv.submat(0, 0, nx - 1, nx - 1) = Varx_inv + Alpha.t() * VmAlpha + Delta.t() * VyDelta;
    // Block 1,2
    ImpCov_inv.submat(0, nx, nx - 1, nx + nm - 1) = -VmAlpha.t() + Delta.t() * VyBeta;
    // Block 1,3
    ImpCov_inv.submat(0, nx + nm, nx - 1, nx + nm + nt - 1) = -VyDelta.t();


    // Block 2,1 = Block 1,2
    ImpCov_inv.submat(nx, 0, nx + nm - 1, nx - 1) = ImpCov_inv.submat(0, nx, nx - 1, nx + nm - 1).t();
    // Block 2,2
    ImpCov_inv.submat(nx, nx, nx + nm - 1, nx + nm - 1) = Varm_inv + Beta.t() * VyBeta;
    // Block 2,3
    ImpCov_inv.submat(nx, nx + nm, nx + nm - 1, nx + nm + nt - 1) = -VyBeta.t();

    // Block 3,1 = Block 1,3
    ImpCov_inv.submat(nx + nm, 0, nx + nm + nt - 1, nx - 1) = ImpCov_inv.submat(0, nx + nm, nx - 1, nx + nm + nt - 1).t();
    // Block 3,2 = Block 2,3
    ImpCov_inv.submat(nx + nm, nx, nx + nm + nt - 1, nx + nm - 1) = ImpCov_inv.submat(nx, nx + nm, nx + nm - 1, nx + nm + nt - 1).t();
    // Block 3,3
    ImpCov_inv.submat(nx + nm, nx + nm, nx + nm + nt - 1, nx + nm + nt - 1) = Vary_inv;

    return ImpCov_inv;
}

double grad_Amat(int index_i,
        int index_j,
        arma::mat & Alpha,
        arma::mat & Beta,
        arma::mat & Delta,
        arma::mat & SampCov,
        arma::mat & S,
        arma::mat & Varx_inv, 
        arma::mat & Varm_inv,
        arma::mat & Vary_inv){

    // compute gradient of loss function for theta(i,j),
    // where theta(i,j) is parameter in asymmetric matrix A

    double grad = 0.0;

    arma::mat B = compute_B(Alpha, Beta, Delta);

    arma::mat ImpCov_inv = inverse_ImpCov(Alpha, Beta, Delta, 
                                             Varx_inv, Varm_inv, Vary_inv);
    arma::mat Ident = eye(size(SampCov));

    arma::mat C = Ident - ImpCov_inv * SampCov;
    
    // note that E is implied covar matrix
    arma::mat E = B * S * B.t();
 
    arma::mat deriv15 = col_outerprod_row(B, E, index_i, index_j);
    deriv15 = deriv15 + deriv15.t();
    grad = trace(ImpCov_inv * deriv15 * C);

    return grad;
}

double grad_AmatV2(int index_i,
        int index_j,
        arma::mat & Alpha,
        arma::mat & Beta,
        arma::mat & Delta,
        arma::mat & SampCov,
        arma::mat & Varx_inv,
        arma::mat & Varm_inv,
        arma::mat & Vary_inv,
        arma::mat & Varx,
        arma::mat & Varm,
        arma::mat & Vary){

    // compute gradient of loss function for theta(i,j),
    // where theta(i,j) is parameter in asymmetric matrix A

    double grad = 0.0;

    arma::mat B = compute_B(Alpha, Beta, Delta);

    arma::mat ImpCov_inv = inverse_ImpCov(Alpha, Beta, Delta, 
                                             Varx_inv, Varm_inv, Vary_inv);
    arma::mat Ident = eye(size(SampCov));

    arma::mat C = Ident - ImpCov_inv * SampCov;
    
    arma::mat E  = compute_E(Alpha, Beta, Delta, Varx, Varm, Vary); 
    
    arma::mat deriv15 = col_outerprod_row(B, E, index_i, index_j);
    deriv15 = deriv15 + deriv15.t();
    grad = trace(ImpCov_inv * deriv15 * C);

    return grad;
}


double grad_Smat(int index_i,
        int index_j,
        arma::mat & Alpha,
        arma::mat & Beta,
        arma::mat & Delta,
        arma::mat & SampCov,
        arma::mat & Varx_inv,
        arma::mat & Varm_inv,
        arma::mat & Vary_inv){

    // compute gradient of loss function for theta(i,j),
    // where theta(i,j) is parameter in symmetric matrix S

    double grad = 0.0;

    arma::mat B = compute_B(Alpha, Beta, Delta);
    
    arma::mat ImpCov_inv = inverse_ImpCov(Alpha, Beta, Delta, 
                                             Varx_inv, Varm_inv, Vary_inv);
    arma::mat Ident = eye(size(SampCov));

    arma::mat C = Ident - ImpCov_inv * SampCov;

    arma::mat Bt = B.t();
    arma::mat deriv15 = col_outerprod_row(B, Bt, index_i, index_j);

    grad = trace(ImpCov_inv * deriv15 * C);

    return grad;
}

double penalty(arma::mat &Alpha,
        arma::mat &Beta,
        arma::mat &Delta,
        double lambda) {
    
  double pen_alpha = pow( (Alpha.n_rows * Alpha.n_cols), 0.4);
  double pen_beta  = pow( (Beta.n_rows * Beta.n_cols),   0.4);
  double pen_delta = pow( (Delta.n_rows * Delta.n_cols), 0.4);
   
    // wt'd sum abs values of values in matrices, with weights
    // sqrt of number of param's
    
    double pen = pen_alpha * lambda * accu(abs(Alpha)) + 
                 pen_beta  * lambda * accu(abs(Beta))   + 
                 pen_delta * lambda * accu(abs(Delta));
    return pen;
}

double update_param(double gradient, double param, double step,
        double lambda) {
    double diff = param - step*gradient;
    double param_new = soft_threshold(diff, step * lambda);
    return param_new;
}


int count_df(arma::mat &A) {
    double eps = 0.001;
    int df = 0;
    for (int i = 0; i < A.n_rows; i++) {
        for (int j = 0; j < A.n_cols; j++) {
            if (fabs(A(i, j)) > eps) {
                df++;
            }
        }
    }
    return (df);
}

int count_df_Vary(arma::mat &V){
    double eps = 0.001;
    int df = 0;
    for (int i = 0; i < V.n_rows; i++) {
        for (int j = i; j < V.n_cols; j++) {
            if (fabs(V(i, j)) > eps) {
                df++;
            }
        }
    }
    return (df);
}


arma::mat gchol(arma::mat matrix) {

/*

A symmetric matrix A can be decomposed as LDL', where L is a lower triangular matrix
with 1's on the diagonal, L' is the transpose of L, and D is diagonal. The inverse of
L is also lower-triangular, with 1's on the diagonal. If all elements of D are positive,
then A must be symmetric positive definite (SPD), and the solution can be reduced
the usual Cholesky decomposition U'U where U is upper triangular and U = sqrt(D) L'.

The main advantage of the generalized form is that it admits of matrices that are not of
full rank: D will contain zeros marking the redundant columns, and the rank of A is the
number of non-zero columns. If all elements of D are zero or positive, then A is a
non-negative definite (NND) matrix. The generalized form also has the (quite minor)
numerical advantage of not requiring square roots during its calculation.
subroutine to do Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
**

subroutine to do Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
**
** arguments are:
**     n         the size of the matrix to be factored
**     **matrix  a ragged array containing an n by n submatrix to be factored
**
**  The factorization is returned in the lower triangle, D occupies the
**    diagonal and the upper triangle is left undisturbed.
**    The lower triangle need not be filled in at the start.
**
**
**  If a column is deemed to be redundant, then that diagonal is set to zero.
**
**   Terry Therneau

DJS: return  L * diag(sqrt(D))

*/

  double EPSILON = 0.000000001;     /* <= EPS is considered a zero */

  double temp;
  int  i,j,k;
  double eps, pivot;
  int rank;
  int n =  matrix.n_rows;

  eps =0;
  for (i=0; i<n; i++) {
    if (matrix(i,i) > eps)  eps = matrix(i,i);
    for (j=(i+1); j<n; j++)  matrix(j,i) = matrix(i,j);
  }
  eps *= EPSILON;

  rank =0;
  for (i=0; i<n; i++) {
    pivot = matrix(i,i);
    if (pivot < eps) matrix(i,i) =0;
    else  {
      rank++;
      for (j=(i+1); j<n; j++) {
	temp = matrix(j,i)/pivot;
	matrix(j,i) = temp;
	matrix(j,j) -= temp*temp*pivot;
	for (k=(j+1); k<n; k++) matrix(k,j) -= temp*matrix(k,i);
      }
    }
  }

  vec diag(n);
  for(i = 0; i < n;i++){
    diag(i) = sqrt(matrix(i,i));
    matrix(i,i) = 1.0;
  }

  for(i=0; i<(n-1); i++){
    for(j=i+1;j<n; j++){
      matrix(i,j) = 0.0;
    }
  }

  matrix = matrix * diagmat(diag);

  return  matrix;
}


arma::mat gchol_inv(arma::mat matrix) {
  // find inversse of lower triangular matrix that
  // is a generalized cholesky  L * diag(sqrt(D))

  int n = matrix.n_rows;
  arma::mat inv = matrix;

  int i,j,k;
  double diag;

  for (k=0; k<n; k++){

    if (inv(k,k) > 0.0) {

      diag = inv(k,k);
      for (i=0; i < k; i++){
	inv(k,i) = inv(k,i)/diag;
      }
      for(i=k; i< n; i++){
	inv(i,k) = inv(i,k)/diag;
      }

      for(i=0; i < n; i++){
	if(i == k) continue;
	for(j=0; j < n; j++){
	  if(j == k) continue;
	  inv(i,j) = inv(i,j) - inv(i,k) * inv(k,j) * diag;
	}
      }
      inv(k,k) = - 1.0/diag;
    }
  }
  for(i=0; i < n; i++){
    for(j=0; j <=i; j++){
      inv(i,j) = - inv(i,j);
    }
  }

  return(inv);
}
void print_vec(arma::vec x){

  for(int i = 0; i < x.size(); i++){
    Rcout << x(i)  << ", ";
  }
  return;
}
void print_mat(arma::mat &x){

  for(int i = 0; i < x.n_rows; i++){
    for(int j = 0; j < x.n_cols; j++){
      Rcout <<  x(i,j)  << ", ";
    }
    Rcout << endl;
  }
  Rcout << endl;

  return;
}


//[[Rcpp::export]]

List rcpp_mvregmed(
            arma::mat alpha,
            arma::mat beta,
            arma::mat delta,
            arma::mat varx,
            arma::mat varm,
            arma::mat vary,
            arma::mat sampcov,
            double sample_size,
            double lambda,
            int max_iter,
            int max_iter_inner,
            double tol = 1e-5,
            double vary_step_size = .05,
            double step_multiplier = .5,
            bool verbose = false) {

    // fit structural equation model with lasso penalties on
    // alpha, beta, delta, but no penalties on vary

     //------------ define parameters/variables
     
    // pen below are used to weight lamdba for different number of
    // params for alpha/beta/delta, because without these weights,
    // find alpha's are selected more frequently than betas and deltas
    // presumably because there can be many more alphas than other params
    
  double pen_alpha = pow( (alpha.n_rows * alpha.n_cols), 0.4);
  double pen_beta  = pow( (beta.n_rows * beta.n_cols),   0.4);
  double pen_delta = pow( (delta.n_rows * delta.n_cols), 0.4);
   
    int iter_inner = 0;
    int index_i = 0;
    int index_j = 0;

    double gradient;
     
    double diff2 = 0.0;
    double gdiff = 0.0;
    double grad_delta = 0.0;
    double grad_old = 0.0;
    
    double param_old = 0.0;
    double param_new = 0.0;
    double param_diff = 0.0;
    double loss_majorize = 0.0;
    double loss_old = 0.0;
    double loss_new = 0.0;
    double pen_loss_begin = 0.0;
    double pen_loss_new = 0.0;
    double pen_loss_old = 0.0;
    double pen_loss_end = 0.0;
    double step = 0.5;             // init step size for alpha, beta, delta
    double step_mult_vary = 0.25;   // init step size for vary
   
    bool converge = false;          // logical for global converge
    bool converge_inner = false;    // logical for inner-loop converge

    bool debug1 = false;
    bool debug2 = false;
     
    
    
    //------------- Initial computations before loops
  
    arma::mat ImpCov(size(sampcov), fill::zeros);

    //  S matrix arranged as (varx, 0..0,     0
    //                        0,    varm,     0
    //                        0,    0..0,  vary)

    // S has dim nx + mn + nt
    // nx = number of instrumental variables (exposures, SNPs, other aliases)
    // nm = number of mediators
    // nt = number of traits

    int nx = varx.n_rows;
    int nm = varm.n_rows;
    int nt = vary.n_rows;

    arma::mat S = blockDiag(varx, varm, vary);

  
    arma::mat varx_inv = pinv(varx);
    arma::mat varm_inv = pinv(varm);
    arma::mat vary_inv = pinv(vary);
    
    // B = inv(I-A)
    arma::mat B = compute_B(alpha, beta, delta);

    // compute implied cov
    ImpCov = B * S * B.t();

    double logdet_varx = logdet_Var(varx);
    double logdet_varm = logdet_Var(varm);
    double logdet_sum = logdet_varx + logdet_varm;
    double logdet_ImpCov = logdet_sum + logdet_Var(vary);

    arma::mat ImpCov_inv = inverse_ImpCov(alpha, beta, delta, varx_inv, varm_inv, vary_inv);
    loss_old = logdet_ImpCov + trace(sampcov * ImpCov_inv);

    pen_loss_old = loss_old + penalty(alpha, beta, delta, lambda);
    pen_loss_begin = pen_loss_old;

    //----------------  outer loop for all parameters
    
    // init for outer loop
    
    int iter = 0;         
    converge = false;

    while ((iter <= max_iter) & !converge) {
            
        iter ++;
        if(debug1) {
            Rcout << "iter = " << iter << endl;
            Rcout << "alpha:" << endl;
            print_mat(alpha);
        }
        
       
        if (verbose & ((iter % 100) == 0)) Rcout << "iter = " << iter << endl;

        //============================== update alpha's ======================================//

        for (int i_alpha = 0; i_alpha < alpha.n_rows; i_alpha++) {
            for (int j_alpha = 0; j_alpha < alpha.n_cols; j_alpha++) {

                if (debug1) {
                    Rcout << "alpha inner loop" << endl;
                    Rcout << "======= alpha ===========" << endl;
                    print_mat(alpha);
                    //Rcout << "======= beta ===========" << endl;
                    //print_mat(beta);
                    //Rcout << "======= delta ===========" << endl;
                    //print_mat(delta);
                    //Rcout << "======= vary ===========" << endl;
                   //print_mat(vary);
                    Rcout << "begin inner loop loss old = " << loss_old << endl;
                }
                
                // inner loop to optimize over a single alpha parameter

                iter_inner = 0;
                converge_inner = false;
                while ((iter_inner < max_iter_inner) &&  !converge_inner) {

                    param_old = alpha(i_alpha, j_alpha);
                    
                    
                    // index_i, index_j are indices where alpha's occur in the
                    // Asymetric A matrix; nx = number of "exposure" variables
                    
                    index_i = nx + i_alpha;
                    index_j = j_alpha;

		    // update gradient of A matrix wrt alpha param
 
                    // old gradient =  grad_Amat(index_i, index_j, alpha, beta, delta, sampcov, S, varx_inv, varm_inv, vary_inv);
                    gradient= grad_AmatV2(index_i, index_j, alpha, beta, delta, sampcov, varx_inv, varm_inv, vary_inv,
                                             varx, varm, vary);
                          
                    if (debug1) {
                        Rcout << "alpha[" << i_alpha << ", " << j_alpha << "] old: ";
                        Rcout << " param " << param_old << ", loss = " << loss_old << ", grad = " << gradient << endl;
                    }
        
                    // optimize step size
                    
                    step = 0.5;
                    for (int i = 0; i < 10; i++) {

                        param_new = update_param(gradient, param_old, step, lambda*pen_alpha);
                        alpha(i_alpha, j_alpha) = param_new;

                        // When alpha is updated, the following should be updated:
                        ImpCov_inv = inverse_ImpCov(alpha, beta, delta, varx_inv, varm_inv, vary_inv);
                        loss_new = logdet_ImpCov + trace(sampcov * ImpCov_inv);

                        param_diff = param_new - param_old;
                        gdiff = gradient * param_diff;
                        diff2 = param_diff * param_diff;
                        loss_majorize = loss_old + gdiff + diff2 / (2.0 * step);

                        if (debug1) {
                            Rcout << "step = " << step << ", alpha[" << i_alpha << ", " << j_alpha << "]";
                            Rcout << " = " << param_new << ", loss new = " << loss_new;
                            Rcout << ", loss maj = " << loss_majorize;
                            Rcout << ", diff^2/(2*step) = " << diff2 / (2.0 * step)  << endl;
                        }

                        if (loss_new <= loss_majorize) {
                            if(debug1) Rcout << "alpha break" << endl;
                            break;
                        }

                        step = step_multiplier * step;

                    } // end of step optimization

                    // if after step opt, loss_new > loss_majorize revert to parma_old
                    if(loss_new > loss_majorize){
                        alpha(i_alpha, j_alpha) = param_old;
                    }
                

                    // When alpha is updated, the following should be updated:
                    ImpCov_inv = inverse_ImpCov(alpha, beta, delta, varx_inv, varm_inv, vary_inv);

                    loss_new = logdet_ImpCov + trace(sampcov * ImpCov_inv);

                    // check convergence
                    pen_loss_new = loss_new + penalty(alpha, beta, delta, lambda);

                    if (fabs(pen_loss_new - pen_loss_old) < tol * (fabs(pen_loss_old) + 1.0)) {
                        converge_inner = true;
                    }

                    pen_loss_old = pen_loss_new;
                    loss_old = loss_new;
                    iter_inner++;

                } // end inner loop

            }
        }

    
        //============================== update beta's ======================================//

        for (int i_beta = 0; i_beta < beta.n_rows; i_beta++) {
            for (int j_beta = 0; j_beta < beta.n_cols; j_beta++) {

                if (debug2) {
                    Rcout << "beta inner loop" << endl;
                    Rcout << "======= alpha ===========" << endl;
                    print_mat(alpha);
                    Rcout << "======= beta ===========" << endl;
                    print_mat(beta);
                    Rcout << "======= delta ===========" << endl;
                    print_mat(delta);
                    Rcout << "======= vary ===========" << endl;
                    print_mat(vary);
                    Rcout << "begin inner loop loss old = " << loss_old << endl;
                }
                
                // inner loop to optimize over a single beta parameter

                iter_inner = 0;
                converge_inner = false;
                while ((iter_inner < max_iter_inner) &&  !converge_inner) {

                    param_old = beta(i_beta, j_beta);
                    // index_i, index_j indices for beta's in Asymmetric A matrix
                    index_i = nx + nm + i_beta;
                    index_j = nx + j_beta;

                    // old gradient = grad_Amat(index_i, index_j, alpha, beta, delta, sampcov, varx_inv, varm_inv, vary_inv);
                    gradient= grad_AmatV2(index_i, index_j, alpha, beta, delta, sampcov, varx_inv, varm_inv, vary_inv,
                                             varx, varm, vary);
                  
                    if (debug2) {
                        Rcout << "beta[" << i_beta << ", " << j_beta << "] old: ";
                        Rcout << " = " << param_old << ", loss = " << loss_old << ", grad = " << gradient << endl;
                    }
                    
                    // optimize step size
                    step = 0.5;
                    for (int i = 0; i < 10; i++) {

                        param_new = update_param(gradient, param_old, step, lambda*pen_beta);
                        beta(i_beta, j_beta) = param_new;

                        // When beta is updated, the following should be updated:
                        ImpCov_inv = inverse_ImpCov(alpha, beta, delta, varx_inv, varm_inv, vary_inv);
                        loss_new = logdet_ImpCov + trace(sampcov * ImpCov_inv);

                        param_diff = param_new - param_old;
                        gdiff = gradient * param_diff;
                        diff2 = param_diff * param_diff;
                        loss_majorize = loss_old + gdiff + diff2 / (2.0 * step);

                        if (debug2) {
                            Rcout << "step beta[" << i_beta << ", " << j_beta << "] ";
                            Rcout << " = " << param_new << ", loss new = " << loss_new << ", loss maj = " << loss_majorize << endl;
                        }
                        if (loss_new <= loss_majorize) {
                            if(debug2) Rcout << "beta break" << endl;
                            break;
                        }

                        step = step_multiplier * step;

                    } // end of step optimization

                    
                    // if after step opt, loss_new > loss_majorize revert to parma_old
                    if(loss_new > loss_majorize){
                        beta(i_beta, j_beta) = param_old;
                    }
     
                    // When beta is updated, the following should be updated:
                    ImpCov_inv = inverse_ImpCov(alpha, beta, delta, varx_inv, varm_inv, vary_inv);
                    loss_new = logdet_ImpCov + trace(sampcov * ImpCov_inv);

                    // check convergence
                    pen_loss_new = loss_new + penalty(alpha, beta, delta, lambda);

                    if (fabs(pen_loss_new - pen_loss_old) < tol * (fabs(pen_loss_old) + 1.0)) {
                        converge_inner = true;
                    }

                    pen_loss_old = pen_loss_new;
                    loss_old = loss_new;
                    iter_inner++;

                } // end inner loop

            }
        }
   
   
        //============================== update delta's ======================================//

        for (int i_delta = 0; i_delta < delta.n_rows; i_delta++) {
            for (int j_delta = 0; j_delta < delta.n_cols; j_delta++) {

                if (debug2) {
                    Rcout << "delta inner loop" << endl;
                    Rcout << "======= alpha ===========" << endl;
                    print_mat(alpha);
                    Rcout << "======= beta ===========" << endl;
                    print_mat(beta);
                    Rcout << "======= delta ===========" << endl;
                    print_mat(delta);
                    Rcout << "======= vary ===========" << endl;
                    print_mat(vary);
                    Rcout << "begin inner loop loss old = " << loss_old << endl;
                }
                
                iter_inner = 0;
                converge_inner = false;
                while ((iter_inner < max_iter_inner) &&  !converge_inner) {

                    param_old = delta(i_delta, j_delta);
                    
                    // index_i, index_j are indices for where delta's occur
                    // in Asymmetric A matrix
                    
                    index_i = nx + nm + i_delta;
                    index_j = j_delta;
                    
                    // gradient of A matrix wrt to delta's
                    // old gradient = grad_Amat(index_i, index_j, alpha, beta, delta, sampcov, varx_inv, varm_inv, vary_inv);
                    gradient= grad_AmatV2(index_i, index_j, alpha, beta, delta, sampcov, varx_inv, varm_inv, vary_inv,
                                             varx, varm, vary);
                  
                    if (debug2) {
                        Rcout << "delta[" << i_delta << ", " << j_delta << "] old: ";
                        Rcout << " param " << param_old << ", loss = " << loss_old << ", grad = " << gradient << endl;
                    }
                    
                    // optimize step size

                    step = 0.5;
                    for (int i = 0; i < 10; i++) {

                        param_new = update_param(gradient, param_old, step, lambda*pen_delta);
                        delta(i_delta, j_delta) = param_new;

                        // When delta is updated, the following should be updated:
                        ImpCov_inv = inverse_ImpCov(alpha, beta, delta, varx_inv, varm_inv, vary_inv);
                        loss_new = logdet_ImpCov + trace(sampcov * ImpCov_inv);

                        param_diff = param_new - param_old;
                        gdiff = gradient * param_diff;
                        diff2 = param_diff * param_diff;
                        loss_majorize = loss_old + gdiff + diff2 / (2.0 * step);

                        if (debug2) {
                            Rcout << "step delta[" << i_delta << ", " << j_delta << "] ";
                            Rcout << " = " << param_new << ", loss new = " << loss_new << ", loss maj = " << loss_majorize << endl;
                        }
                        
                        if (loss_new <= loss_majorize) {
                            if(debug2) Rcout << "delta break" << endl;
                            break;
                        }

                        step = step_multiplier * step;

                    } // end of step optimization

                    // if after step opt, loss_new > loss_majorize revert to parma_old
                    if(loss_new > loss_majorize){
                        delta(i_delta, j_delta) = param_old;
                    }
     
                    // When delta is updated, the following should be updated:
                   ImpCov_inv = inverse_ImpCov(alpha, beta, delta, varx_inv, varm_inv, vary_inv);
                   loss_new = logdet_ImpCov + trace(sampcov * ImpCov_inv);

                    // check convergence
                    pen_loss_new = loss_new + penalty(alpha, beta, delta, lambda);

                    if (fabs(pen_loss_new - pen_loss_old) < tol * (fabs(pen_loss_old) + 1.0)) {
                        converge_inner = true;
                    }

                    pen_loss_old = pen_loss_new;
                    loss_old = loss_new;
                    iter_inner++;

                } // end inner loop for delta
            }
        }

    
    //============================== update vary's ======================================//
      
        // note that loops cover upper triangle of vary

        for (int i_vary = 0; i_vary < vary.n_rows; i_vary++) {
            for (int j_vary = i_vary; j_vary < vary.n_cols; j_vary++) {

                if (debug2) {
                    Rcout << "vary inner loop" << endl;
                    Rcout << "======= alpha ===========" << endl;
                    print_mat(alpha);
                    Rcout << "======= beta ===========" << endl;
                    print_mat(beta);
                    Rcout << "======= delta ===========" << endl;
                    print_mat(delta);
                    Rcout << "======= vary ===========" << endl;
                    print_mat(vary);
                    Rcout << "begin inner loop loss old = " << loss_old << endl;
                }
                
                // inner loop to optimize over a single parameter

                iter_inner = 0;
                converge_inner = false;
                
                while ((iter_inner < max_iter_inner) &&  !converge_inner) {
                    
                    iter_inner++;
                    
                    param_old = vary(i_vary, j_vary);
                    
                    // index_i, index_j are indices for where vary's occur
                    // in symmetric S matrix
                    
                    index_i = nx + nm + i_vary;
                    index_j = nx + nm + j_vary;
                    
                    // gradient of S mat wrt to vary
                    gradient = grad_Smat(index_i, index_j, alpha, beta, delta, sampcov, varx_inv, varm_inv, vary_inv);

                    grad_delta = fabs(gradient - grad_old)/ (1.0 + fabs(grad_old));
                    
		    // if large change in gradient, revert to old parm and jump out of iter_inner loop
                    // is following portable?
                    if(isNaN(gradient)){
                        break;
                    }
		    if( fabs(grad_delta) > 10.0){
                        break;
		    }
                    grad_old = gradient;
               
                    if (debug2) {
                        Rcout << "var[" << i_vary << ", " << j_vary << "] old: ";
                        Rcout << " param " << param_old << ", loss = " << loss_old << ", grad = " << gradient << endl;
                    }
                    
                    // optimize step size
                    step = vary_step_size;
                    for (int i = 0; i <10; i++) {
                        
                        // note no penalty on vary
                        param_new = param_old - step*gradient;
                     
                        // bound diag var away from 0
                        if ((i_vary == j_vary) && param_new < 0.00001) {
                            step = step_mult_vary * step;
                            continue;
                        }
                                
                        // an item of vary is updated, so update anything depending
                        // on vary

                        vary(i_vary, j_vary) = param_new;
                        vary(j_vary, i_vary) = param_new;
			
		
			//Rcout << "before pinv(vary) in mvregmed - 1" << endl;
			//Rcout << "grad = " << gradient << ", step = " << step <<  endl;
			//print_mat(vary);

                        vary_inv = pinv(vary);
                         
                        logdet_ImpCov = logdet_sum + logdet_Var(vary);

                        ImpCov_inv = inverse_ImpCov(alpha, beta, delta, varx_inv, varm_inv, vary_inv);
                        loss_new = logdet_ImpCov + trace(sampcov * ImpCov_inv);
                        
                        if (debug2) {
                            Rcout << "step  var[" << i_vary << ", " << j_vary << "] ";
                            Rcout << " = " << param_new << ", loss = " << loss_new << endl;
                            Rcout << "step logdet_sum = " << logdet_sum << ", logdet(vary) = " << logdet_Var(vary);
                            Rcout << ", tr = " << trace(sampcov * ImpCov_inv) << endl;
                        }
                        
                        if ( loss_new <= loss_old) {
                            if(debug2) Rcout << "vary break" << endl;  
                            break;
                        }
                        step = step_mult_vary * step;
                    }// end of step optimization
                   
                   
                    // if after step opt, loss_new > loss_old, revert to loss_old and parm_old
                    if(loss_old < loss_new){
                        vary(i_vary, j_vary) = param_old;
                        vary(j_vary, i_vary) = param_old;
			//Rcout << "before pinv(vary) in mvregemd - 2" << endl;
                        vary_inv = pinv(vary);
                        logdet_ImpCov = logdet_sum + logdet_Var(vary);
                        ImpCov_inv = inverse_ImpCov(alpha, beta, delta, varx_inv, varm_inv, vary_inv);
                        loss_new = logdet_ImpCov + trace(sampcov * ImpCov_inv);    
                    }
             

                    // check convergence
                    pen_loss_new = loss_new + penalty(alpha, beta, delta, lambda);
 
                    if (fabs(pen_loss_new - pen_loss_old) < tol * (fabs(pen_loss_old) + 1.0)) {
                        converge_inner = true;
                    }

                    pen_loss_old = pen_loss_new;
                    loss_old = loss_new;
                   // moved to top of inner loop  iter_inner++;
             
                } // end inner loop

            }
        }

        // this completes sequentially optimizing over each parameter
        // check for global convergence by comparing the penalized loss function
        // at the beginning of a global iter to the end of a global iter

        if (fabs(pen_loss_new - pen_loss_begin) < tol * (fabs(pen_loss_new) + 1.0)) {
            converge = true;
        }

        pen_loss_begin = pen_loss_new;
        
    } 
    // end outer loop

    // compute BIC

    int df_alpha = count_df(alpha);
    int df_beta  = count_df(beta);
    int df_delta = count_df(delta);
    int df_vary  = count_df_Vary(vary);
    int df = df_alpha + df_beta + df_delta + df_vary;
    
    double bic = loss_new * sample_size + log(sample_size) * double(df);

    if(debug1){
        Rcout << "alpha matrix before return:" << endl;
        print_mat(alpha);
    }

    return Rcpp::List::create(Rcpp::Named("alpha") = alpha,
                              Rcpp::Named("beta") = beta,
                              Rcpp::Named("delta") = delta,
                              Rcpp::Named("varx") = varx,
                              Rcpp::Named("varm") = varm,
                              Rcpp::Named("vary") = vary,                         
                              Rcpp::Named("lambda") = lambda,
                              Rcpp::Named("converge") = converge,
                              Rcpp::Named("iter") = iter,
                              Rcpp::Named("loss") = loss_new,
                              Rcpp::Named("penloss") = pen_loss_new,
                              Rcpp::Named("bic") = bic,
                              Rcpp::Named("df") = df,
                              Rcpp::Named("df.alpha") = df_alpha, 
                              Rcpp::Named("df.beta") = df_beta,
                              Rcpp::Named("df.delta") = df_delta,
                              Rcpp::Named("df.vary") = df_vary);
}

