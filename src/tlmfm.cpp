// This script implement linear matrix factor model based on variational inference.
// Date: 2024-07-07

// Revised log:


#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include<ctime>
//#include<boost/math/tools/minima.hpp>

// #define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;
// using boost::math::tools::brent_find_minima;


// Define global variables
//double X_count_ij, a_i, invLambda_j,Mu_x_ij;
/*
 * Auxiliary
 */
//' @keywords internal
 //' @noRd
 //'
 // diag(W0* Cki * W0^t)
 vec decomp(const mat& Cki, const mat& W0){
   vec s, tmp1;
   mat U, V, WC12;
   svd(U, s, V, Cki);
   WC12 = W0 * (U * diagmat(sqrt(s))); // p * q
   tmp1 = sum(WC12 % WC12, 1);
   return tmp1;
 }
List irlbaCpp(const mat& X, const int& q){
  Rcpp::Environment irlba("package:irlba");
  Rcpp::Function f = irlba["irlba"];
  
  return f(Named("A") = X, Named("nv") = q);
  
}  

double afun(const mat& Yt, const mat& mu, const mat& R, const mat& C, 
            const vec& Lam1, const vec& Lam2, const double& nu,
            const mat& Mt, const mat& S1t, const mat& S2t){
  int r1=S1t.n_rows, r2 = S2t.n_rows, p1 = Yt.n_rows, p2=Yt.n_cols;
  double a = 0.0;
  mat tmp_mat;
  mat CtC = C.t()* (C % repmat(1.0/Lam2, 1, r2));
  mat RtR = R.t()*(R % repmat(1.0/Lam1, 1, r1));
  mat S1oS2(r1*r2, r1*r2, fill::zeros); // need to be speeded up!!!!!!!!!
  tmp_mat = (Yt - mu - R*Mt*C.t()) % repmat(1.0/sqrt(Lam1), 1, p2) % repmat(1.0/sqrt(Lam2.t()), p1, 1);
  a += accu(tmp_mat % tmp_mat)/nu;
  // logPy += trace(diagmat(1.0/Lam1)* (Y.slice(t) - mu - R*M.slice(t)*C.t()) * diagmat(1.0/Lam2) * (Y.slice(t) - mu - R*M.slice(t)*C.t()).t());
  a += trace(CtC * S2t) * trace(RtR*S1t) / nu + 1; 
  
  return a;
}

vec atvec_fun(const cube& Y, const cube& M, const mat& mu, const cube& S1,
                       const cube& S2, const mat& R, const mat& C,
                       const vec& Lam1, const vec& Lam2,const double& nu){
  int t, T = Y.n_slices;
  vec atvec(T);
  for(t=0; t<T; ++t){
    atvec(t) = afun( Y.slice(t), mu, R, C, Lam1, Lam2, nu, M.slice(t), S1.slice(t), S2.slice(t));
  }
  
  return atvec;
}

void update_mu(const cube& Y, const cube& M, const mat& R, const mat& C,
                 const vec& atvec, mat& mu){
  
  int t, r1 = R.n_cols, r2 = C.n_cols, p1=mu.n_rows, p2=mu.n_cols, T = Y.n_slices;
  double sum_at = accu(1.0/atvec);
  mat mY(p1,p2, fill::zeros), mM(r1, r2, fill::zeros);
  for(t=0; t<T; ++t){
    mY += Y.slice(t) /atvec(t);
    mM += M.slice(t) /atvec(t);
  }
  mu = (mY - R * mM * C.t())/ sum_at;
}


void update_R(const cube& Y, const cube& M, const mat& mu, const cube& S1,
              const cube& S2, const mat& C,
              const vec& Lam2, const vec& atvec, mat& R){
  int t, r1 = S1.n_rows, r2 = S2.n_rows, p1=mu.n_rows, T = Y.n_slices;
  mat CLCt = (C % repmat(1/Lam2, 1, r2)).t() * C; // r2*r2
  mat sMt = zeros(r1, r1), sSt1 = zeros(r1, r1); // can be speeded up using trace
  mat sHLCM = zeros(p1, r1);
  double at=0;
  for(t=0; t<T; ++t){
    at = atvec(t);
    sMt += M.slice(t) * CLCt * M.slice(t).t() / at;
    sSt1 += trace(S2.slice(t)* CLCt) * S1.slice(t)/ at; // r1*r1
    sHLCM += (Y.slice(t)- mu)* (repmat(1/Lam2, 1, r2) % C) * M.slice(t).t() / at; // p1* r1
  }
  // Rprintf("Good update for R1!\n");
  R = sHLCM * (sMt+ sSt1).i();
}


void update_C(const cube& Y, const cube& M, const mat& mu, const cube& S1,
              const cube& S2, const mat& R,
              const vec& Lam1, const vec& atvec, mat& C){
  int t, r1 = S1.n_rows, r2 = S2.n_rows, p2=mu.n_cols, T = Y.n_slices;
  mat RLRt = (R% repmat(1/Lam1, 1, r1)).t() * R; // r1*r1
  mat sMb = zeros(r2, r2), sSt2 = zeros(r2, r2);
  mat sHLRM = zeros(p2, r2);
  double at=0;
  for(t=0; t<T; ++t){
    at = atvec(t);
    sMb += M.slice(t).t() * RLRt * M.slice(t) / at; // r2*r2
    sSt2 += trace(S1.slice(t)*RLRt) * S2.slice(t) / at; // r1*r1
    // M.slice(t).print();
    sHLRM += (Y.slice(t)- mu).t()* (repmat(1/Lam1, 1, r1) % R) * M.slice(t) / at; // p2* r2
  }

  C = sHLRM * (sMb+ sSt2).i();
}



void update_Lam1(const cube& Y,  const cube& M, const mat& mu,
                 const mat& R, const mat& C, const cube& S1,
                 const cube& S2, const vec& Lam2, const double& nu,
                 const vec& atvec, vec& Lam1){
  int t, p1=mu.n_rows, p2 = mu.n_cols, r2 = C.n_cols, T = Y.n_slices;
  int p = p1*p2;
  vec Lam_tmp(p1, fill::zeros), vec_tmp;
  mat CtC = C.t()* (C % repmat(1.0/Lam2, 1, r2)); // C \in p2*r2
  mat mat_tmp;
  double at=0.0;
  for(t=0; t<T; ++t){
    at = atvec(t);
    mat_tmp = (Y.slice(t) - mu - R*M.slice(t)*C.t()); // %;
    vec_tmp = decomp(S1.slice(t), R);
    Lam_tmp += vec_tmp * trace(CtC * S2.slice(t)) / (T*p2*at*nu); // p1 vector
    Lam_tmp += mean(mat_tmp %  repmat(1.0/Lam2.t(), p1, 1) % mat_tmp, 1) /(T*at*nu);
  }
  
  Lam1 = Lam_tmp*(nu+p);
}

void update_Lam2(const cube& Y,  const cube& M, const mat& mu,
                 const mat& R,const mat& C, const cube& S1,
                 const cube& S2,const vec& Lam1,  const double& nu,
                 const vec& atvec, vec& Lam2){
  int t, p1=mu.n_rows, p2 = mu.n_cols, r1 = R.n_cols, T = Y.n_slices;
  int p = p1*p2;
  rowvec Lam_tmp(p2, fill::zeros);
  mat RtR = R.t() * (R % repmat(1.0/Lam1, 1, r1));
  mat mat_tmp;
  vec vec_tmp;
  double at=0.0;
  for(t=0; t<T; ++t){
    at = atvec(t);
    mat_tmp = (Y.slice(t) - mu - R*M.slice(t)*C.t()); // %;
    vec_tmp = decomp(S2.slice(t), C);
    Lam_tmp += vec_tmp.t() * trace(RtR * S1.slice(t)) / (T*p1*at*nu); // p2 rowvector
    Lam_tmp += mean(mat_tmp %  repmat(1.0/Lam1, 1, p2) % mat_tmp, 0)/(T*at*nu) ;
  }
  Lam2 = Lam_tmp.t()*(nu+p);
}
double logCp(const double& nu, const int& p){
  
  double p1 = (p+nu)/2.0;
  // Rprintf("Gamma(p+nu/2)= %4f \n", lgamma(p1));
  // Rprintf("Gamma(nu/2)= %4f \n", lgamma(nu/2));
  double y = -p/2*log(nu) + lgamma(p1) - lgamma(nu/2); // log(Gamma(x)) in std library.
  return y;
}

double objfun_nu(const cube& Y,  const mat& mu, const mat& R,
                 const mat& C,  const vec& Lam1, const vec& Lam2, const double& nu,
                 const cube& M, const cube& S1, const cube& S2){
  int t,  p1 = Y.n_rows, p2=Y.n_cols,  T = Y.n_slices;
  int p = p1*p2;
  double obj = T* logCp(nu, p);
  double tmp_val = 0.0;
  for(t=0; t<T; ++t){
    tmp_val = afun( Y.slice(t), mu, R, C, Lam1, Lam2, nu, M.slice(t), S1.slice(t), S2.slice(t));
    obj -= (nu+p)/2 * log(tmp_val);
  }
  return obj;
}

void update_nu(const cube& Y,  const mat& mu, const mat& R,
               const mat& C,  const vec& Lam1, const vec& Lam2,
               const cube& M, const cube& S1, const cube& S2,
               const vec& nu_grid, double& nu){
  int i, n_grid = nu_grid.n_elem;
  vec obj_vec(n_grid);
  for(i=0; i<n_grid; ++i){
    obj_vec(i) = objfun_nu( Y, mu, R, C, Lam1, Lam2, nu_grid(i), M, S1, S2);
  }
  // Rprintf("obj_vec= %4f \n");
  // obj_vec.t().print();
  nu= nu_grid(index_max(obj_vec));
  
}

void add_IC_Orth(mat& B){
  // Try the orthogonal matrix method
  int qs1 = B.n_cols, p = B.n_rows;
  mat U1, V1;
  vec s1;
  svd(U1, s1, V1, B);
  vec signU1 = sign(U1.row(0).t());
  B = sqrt(p) * U1.cols(0,qs1-1) * diagmat(signU1.subvec(0,qs1-1));
}

void reupdate_F(const cube& Y, const mat& mu, const mat& R,
           const mat& C, cube& F){
  int t, p1 = Y.n_rows, p2=Y.n_cols,  T = Y.n_slices;
  for(t=0; t<T; ++t){
    F.slice(t) = R.t()* (Y.slice(t)-mu) * C /(p1 * p2);
  }
}


double calELBO(const cube& Y,  const mat& mu, const mat& R,
               const mat& C,  const vec& Lam1, const vec& Lam2, const double& nu,
              const cube& M, const cube& S1, const cube& S2){

  int t, r1=S1.n_rows, r2 = S2.n_rows, p1 = Y.n_rows, p2=Y.n_cols,  T = Y.n_slices;
  int p = p1*p2;
  double cp = logCp(nu, p);
  // Rprintf("cp= %4f \n", cp);
  // logP(Y|F)
  double logPy = T* p2*accu(log(Lam1)) + T* p1 *accu(log(Lam2))-2*T*cp;
  double tmp_val = 0.0;
  for(t=0; t<T; ++t){
    tmp_val = afun( Y.slice(t), mu, R, C, Lam1, Lam2, nu, M.slice(t), S1.slice(t), S2.slice(t));
    logPy += (nu+p) * log(tmp_val);
  }
  logPy = -0.5*logPy;

  //log P(F)
  double logPf = 0.0;
  for(t=0; t<T; ++t){
    
    logPf += accu(M.slice(t) % M.slice(t)) + trace(S2.slice(t)) * trace(S1.slice(t));
  }
  logPf = -0.5*logPf;

  // entropy
  double entropy=0.0;
  for(t=0; t<T; ++t){
    entropy += r1*log(det(S2.slice(t))) + r2*log(det(S1.slice(t)));
  }
  entropy = entropy*0.5;

  double ELBO = logPy + logPf + entropy;
  return ELBO;
}



void VB_Estep(const cube& Y, const mat& mu, const mat& R,
              const mat& C, const vec& Lam1, const vec& Lam2, const int& nu,
               cube& M, cube& S1, cube& S2){

  int t, r1=S1.n_rows, r2 = S2.n_rows, p1 = mu.n_rows, p2 = mu.n_cols, T = Y.n_slices;
  int p = p1*p2;
  // Rprintf("Good E-step1!\n");
  // double elbo1 = calELBO(Y, mu, R, C, Lam1, Lam2, nu, M, S1, S2);
  // update M, S1 and S2
  mat tL1 = R.t()* (R % repmat(1.0/Lam1, 1, r1)); // cache
  mat tL2 = C.t() * (C % repmat(1.0/Lam2, 1, r2));
  mat tmpMat_Mt = (nu+p)/(nu)*kron(tL2, tL1);
  mat tmpVec_Mt;
  // Rprintf("Good E-step2!\n");
  double at = 0.0;
  for(t=0; t<T; ++t){
    // update M
    at = afun(Y.slice(t), mu, R, C, Lam1, Lam2, nu, M.slice(t), S1.slice(t), S2.slice(t));
    tmpVec_Mt = (tmpMat_Mt/at + eye(r1*r2, r1*r2)).i() * ((R % repmat(1.0/Lam1, 1, r1)).t()* (Y.slice(t) - mu)  //
                                 * (repmat(1.0/Lam2, 1, r2)% C)/(nu*at)*(nu+p)).as_col();
    tmpVec_Mt.reshape(r1, r2); // only update M
    M.slice(t) = tmpVec_Mt;
    // update S
    S1.slice(t) = (nu+p)/(nu*at)* trace(S2.slice(t) * tL2)* tL1 + trace(S2.slice(t))*eye(r1,r1);
    S1.slice(t) = r2*  S1.slice(t).i();
    S2.slice(t) = (nu+p)/(nu*at)* trace(S1.slice(t) * tL1)* tL2 + trace(S1.slice(t))*eye(r2,r2);
    S2.slice(t) = r1*  S2.slice(t).i();
  }
  // double elbo3 = calELBO(Y, mu, R, C, Lam1, Lam2, nu,  M, S1, S2);
  // Rprintf("dF= %4f \n", elbo3 - elbo1);
}




// [[Rcpp::export]]
Rcpp::List VB_tLMFMcpp(const arma::cube& Y,
                      const arma::cube& M_int, const arma::cube& S1_int,
                      const arma::cube& S2_int, const arma::mat& mu_int,
                      const arma::mat& R_int, const arma::mat& C_int,
                      const arma::vec& Lam1_int, const arma::vec& Lam2_int,
                      const double& nu_int, const bool& add_IC_iter,
                      const double& epsELBO, const int& maxIter,  
                      const bool& verbose,const bool& cal_eigs=false){

  // whether add nvec for Binomial variables???? Now assume n_j = 1 for all.
  // int  T = Y.n_slices;
  vec nu_grid = linspace(1, 60, 30);
  // Rprintf("Initialize!\n");
  // Initialize
  mat mu(mu_int), R(R_int),  C(C_int);
  cube M(M_int), S1(S1_int), S2(S2_int);
  vec Lam1(Lam1_int), Lam2(Lam2_int);
  double nu(nu_int);

  vec ELBO_vec(maxIter), atvec;
  ELBO_vec(0) = -1e20;
  int iter;


  for(iter = 1; iter < maxIter; ++iter){

    atvec = atvec_fun( Y, M, mu, S1, S2, R, C, Lam1, Lam2, nu);
    // Rprintf("E step starting!\n");
    // VB E-step
    VB_Estep(Y, mu, R, C,  Lam1, Lam2, nu, M, S1, S2);
    // Rprintf("Finish E step!\n");
    //VB M-step

    // double elbo1 = calELBO(Y, mu, R, C,  Lam1, Lam2, nu,  M, S1, S2);
    //update Mu
    // Rprintf("update Mu \n");
    update_mu( Y, M, R, C, atvec, mu);
    // double elbo2 = calELBO(Y, mu, R, C,  Lam1, Lam2, nu, M, S1, S2);
    // Rprintf("dmu= %4f \n", elbo2 - elbo1);

    // update R
    // Rprintf("update R \n");
    update_R( Y, M, mu, S1, S2, C, Lam2, atvec, R) ;
    // double elbo3 =  calELBO(Y, mu, R, C,  Lam1, Lam2, nu,  M, S1, S2);
    // Rprintf("dR= %4f\n", elbo3 - elbo2);

    // update C
    // Rprintf("update C\n");
    update_C( Y, M, mu, S1, S2, R, Lam1, atvec, C) ;
    // double elbo4 = calELBO(Y, mu, R, C,  Lam1, Lam2, nu,  M, S1, S2);
    // Rprintf("dC= %4f \n", elbo4 - elbo3);
    if(add_IC_iter){
      add_IC_Orth(R);
      add_IC_Orth(C);
    }
    // update Lam1
    // Rprintf("update Lambda1\n");
    update_Lam1( Y, M, mu, R, C, S1, S2, Lam2, nu,atvec,  Lam1);
    // double elbo5 = calELBO(Y, mu, R, C,  Lam1, Lam2, nu,  M, S1, S2);
    // Rprintf("dLam1= %4f \n", elbo5 - elbo4);
    
    // update Lambda
    // Rprintf("update Lambda2\n");
    update_Lam2( Y, M, mu, R, C, S1, S2, Lam1,nu, atvec, Lam2);
    // Lam2.subvec(0,4).print();
    // double elbo6 =  calELBO(Y, mu, R, C,  Lam1, Lam2, nu,  M, S1, S2);
    // Rprintf("dLam2= %4f\n", elbo6 - elbo5);
    
    // update nu
    // Rprintf("update nu\n");
    // update_nu( Y, mu, R, C, Lam1, Lam2, M, S1, S2, nu_grid, nu);
    // double elbo7 =  calELBO(Y, mu, R, C,  Lam1, Lam2, nu,  M, S1, S2);
    // Rprintf("dnu= %4f\n", elbo7 - elbo6);


    ELBO_vec(iter) =  calELBO(Y, mu, R, C,  Lam1, Lam2, nu,  M, S1, S2);

    if(verbose){
      Rprintf("iter = %d, ELBO= %4f, dELBO=%4f \n",
              iter +1, ELBO_vec(iter), abs(ELBO_vec(iter)  - ELBO_vec(iter-1))/ abs(ELBO_vec(iter-1)));
    }
    if(abs((ELBO_vec(iter)  - ELBO_vec(iter-1))/ ELBO_vec(iter-1)) < epsELBO) break;
  }
  if(!add_IC_iter){
    add_IC_Orth(C);
    add_IC_Orth(R);
    // re-update F after adding identifiability for C and R.
    reupdate_F( Y, mu, R, C, M);
  }
  update_nu( Y, mu, R, C, Lam1, Lam2, M, S1, S2, nu_grid, nu);
  // Select the number of factors
  vec dC, dR;
  int rmax = std::max(R.n_cols, C.n_cols);
  if(cal_eigs){
    //Rprintf("select the number of factors, rmax=%d...\n ", rmax);
    int t, T = Y.n_slices, p1 = R.n_rows, p2 = C.n_rows;
    mat matR(p1, p1, fill::zeros), matC(p2, p2, fill::zeros);
    mat tmpmat1, tmpmat2;
    double at;
    for(t = 0; t<T; ++t){
      at = afun( Y.slice(t), mu, R, C, Lam1, Lam2, nu, M.slice(t), S1.slice(t), S2.slice(t));
      tmpmat1 = (Y.slice(t) - mu)* C;
      matR += tmpmat1*tmpmat1.t()/at;
      tmpmat2 = trans(Y.slice(t) - mu)* R; // weighted covariance.
      matC += tmpmat2*tmpmat2.t()/at;
    }
    Rcpp::List svdC = irlbaCpp(matC, rmax);
    Rcpp::List svdR = irlbaCpp(matR, rmax);
    vec dC1 = svdC["d"];
    vec dR1 = svdR["d"];
    dC = dC1;
    dR = dR1;
  }else{
    vec dC, dR;
    dC = zeros(rmax,1);
    dR = dC;
  }
  // output return value
  List resList = List::create(
    Rcpp::Named("hF") = M,
    Rcpp::Named("hmu") = mu,
    Rcpp::Named("hR") = R,
    Rcpp::Named("hC") = C,
    Rcpp::Named("hnu") = nu,
    Rcpp::Named("hLambda1") = Lam1,
    Rcpp::Named("hLambda2") = Lam2,
    Rcpp::Named("dR") = dR,
    Rcpp::Named("dC") = dC,
    Rcpp::Named("ELBO") = ELBO_vec(iter-1),
    Rcpp::Named("dELBO") = ELBO_vec(iter-1)  - ELBO_vec(iter-2),
    Rcpp::Named("ELBO_seq") = ELBO_vec.subvec(0, iter-1)
  );
  return(resList);

}



