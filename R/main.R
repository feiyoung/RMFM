# generate man files
# devtools::document()
# R CMD check --as-cran RMFM_1.1.0.tar.gz
## usethis::use_data(dat_r2_mac)
# pkgdown::build_site()
# pkgdown::build_home()
# pkgdown::build_reference()
# pkgdown::build_article("COAPsimu")
# pkgdown::build_article("ProFASTdlpfc2")


# Basic functions ---------------------------------------------------------
.logDiffTime <- function(main = "", t1 = NULL, verbose = TRUE, addHeader = FALSE,
                         t2 = Sys.time(), units = "mins", header = "*****",
                         tail = "elapsed.", precision = 3){
  
  # main = ""; t1 = NULL; verbose = TRUE; addHeader = FALSE;
  # t2 = Sys.time(); units = "mins"; header = "###########";
  # tail = "elapsed."; precision = 3
  if (verbose) {
    timeStamp <- tryCatch({
      dt <- abs(round(difftime(t2, t1, units = units),
                      precision))
      if (addHeader) {
        msg <- sprintf("%s\n%s : %s, %s %s %s\n%s",
                       header, Sys.time(), main, dt, units, tail,
                       header)
      }
      else {
        msg <- sprintf("%s : %s, %s %s %s", Sys.time(),
                       main, dt, units, tail)
      }
      if (verbose)
        message(msg)
    }, error = function(x) {
      if (verbose)
        message("Time Error : ", x)
    })
  }
  
  return(invisible(0))
}


.logTime <- function(main='', prefix='*****', versoe=TRUE){
  
  if(versoe){
    message(paste0(Sys.time()," : ", prefix," ",  main))
  }
  
  
}

# Unified functions for tLMFM ----------------------------------------------
mat2orth <- function(mat){
  # mat <- R0
  qr1 <- qr(mat)
  Q <- qr.Q(qr1)
  return(Q * sqrt(nrow(mat)))
}


#' Fit the high-dimensional robust matrix factor model
#' @description Fit the high-dimensional robust matrix factor model via variational inference.
#' @param X a  p1* p2*T array, which is the observed  matrix from each individual, where T is the sample size.
#' @param r1 an optional positive integer, specify the number of row factors; default as 10.
#' @param r2 an optional positive integer, specify the number of column factors; default as 10.
#' @param epsELBO  an optional positive value, tolerance of relative variation rate of the variational lower bound value, default as '1e-9'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed an optional integer, specify the random seed for reproducibility in initialization.
#' @param cal_eigs an optional logical value, specify whether calculate the eigenvalues of covariance matrix, default as \code{FALSE}.
#' @return return a list including the following components:
#' \itemize{
#'   \item \code{hF} - a r1* r2*T array, which is the estimated factor  matrix for each individual, where T is the sample size.
#'   \item \code{hmu} - a p1-by-p2 matrix, the estimated mean matrix.
#'   \item \code{hR} - the estimated row  loading matrix.
#'   \item \code{hC} - the estimated column loading matrix.
#'   \item \code{hnu} - the estimated degree of freedom for the error term.
#'   \item \code{hLambda1} - a p1 vector, the estimated row scatter matrix for error.
#'   \item \code{hLambda2} - a p2 vector, the estimated column scatter matrix for error.
#'   \item \code{dR} - \code{NULL} if \code{cal_eigs=FALSE}; a group of eigenvalues of the  sample covariance across rows if \code{cal_eigs=TRUE}.
#'   \item \code{dC} - \code{NULL} if \code{cal_eigs=FALSE}; a group of eigenvalues of the  sample covariance across columns if \code{cal_eigs=TRUE}.
#'   \item \code{ELBO} -  the ELBO value when algorithm stops;
#'   \item \code{ELBO_seq} - the sequence of ELBO values.
#'   \item \code{time_use} - the running time in model fitting of RMFM;
#' }
#' @details None
#' @seealso None
#' @references None
#' @export
#' @useDynLib RMFM, .registration = TRUE
#' @importFrom  Rcpp evalCpp
#' @importFrom stats  rnorm
#' @examples
#' r1 <- 4; r2 <- 3;
#' Tt <- 100; type <- 'MatrixT'
#' p1 <- 50; p2 <- 50
#' datlist <- gendata_rmfm(i = 1,  Tt = Tt,p1 =p1, p2=p2, r1=r1, r2=r2,
#'                         rho=1, type= 'MatrixT', nu=1)
#' str(datlist)
#' reslist <- RMFM(X=datlist$X, r1=r1, r2=r2,  verbose = TRUE, maxIter = 6)


RMFM <- function(X, r1=10, r2=10,  epsELBO=1e-9, maxIter=30, verbose=TRUE, seed=1,
                  cal_eigs=FALSE){
  
  if(!inherits(X, 'array')){
    stop("RMFM: X must be an three-dimensional array!")
  }
  if(r1<1 || r2<1){
    stop("RMFM: both r1 and r2 must be integters greater than 0!")
  }
  if(!inherits(cal_eigs, 'logical')){
    stop("RMFM: cal_eigs must be a logical value!")
  }
  
  message("Calculate the initial values...")
  Tt <- dim(X)[3]
  p1 <- nrow(X); p2 <- ncol(X)
  
  S1_int <- array(NA, dim=c(r1,r1, Tt+1))
  S2_int <- array(NA, dim=c(r2,r2, Tt+1))
  for(t in 1:(Tt)){
    S1_int[,,t] <- diag(rep(1, r1))
    S2_int[,,t] <- diag(rep(1, r2))
  }
  nu_int <- 5
  Sigma_y_int <- array(1, dim=c(p1,p2, Tt))
  mu_int <- apply(X, c(1,2), mean)
  Lam1_int <- rep(1, p1); Lam2_int <- rep(1, p2)
  
  set.seed(seed)
  R_int <- matrix(rnorm(p1*r1), p1, r1)/(p1); C_int <- matrix(rnorm(p2*r2), p2, r2)/(p2)
  M_int <- array(rnorm(r1*r2*(Tt)), dim=c(r1,r2, Tt))
  message("Finish computing the initial values!")
  tstart <- Sys.time()
  message("Start iterative algorithm...")
  tic <- proc.time()
  reslist <- VB_tLMFMcpp(X, M_int, S1_int, S2_int, mu_int, R_int, C_int, Lam1_int, Lam2_int,
                         nu_int, add_IC_iter = FALSE, epsELBO, maxIter, verbose, cal_eigs) 
  toc <- proc.time()
  .logDiffTime(sprintf(paste0("%s Finish  iterative algorithm"), "*****"), t1 = tstart, verbose = verbose)
  reslist$time.use <- toc[3] - tic[3]
  
  return(reslist)
  
}



#' Select the structure dimension of factor matrix
#' @description Select the structure dimension of factor matrix in the high-dimensional robust matrix factor model
#' @param X a  p1* p2*T array, which is the observed  matrix from each individual, where T is the sample size.
#' @param r_max an optional positive integer, specify the upper bound of row and column factors; default as 10.
#' @param epsELBO  an optional positive value, tolerance of relative variation rate of the variational lower bound value, default as '1e-9'.
#' @param maxIter the maximum iteration of the VEM algorithm. The default is 30.
#' @param verbose a logical value, whether output the information in iteration.
#' @param seed an optional integer, specify the random seed for reproducibility in initialization.
#' @return return a list including the following components:
#' \itemize{
#'   \item \code{rvec} - a two-dimensional vector, the estimated row and column numbers of factors.
#'   \item \code{svrMat} - a r_max-by-2 matrix, the singular value ratios.
#' }
#' @details None
#' @seealso None
#' @references None
#' @importFrom COAP RR_COAP
#' @import irlba
#' @importFrom  irlba irlba
#' @export
#' @examples
#' r1 <- 4; r2 <- 3;
#' Tt <- 100; type <- 'MatrixT'
#' p1 <- 50; p2 <- 50
#' datlist <- gendata_rmfm(i = 1,  Tt = Tt,p1 =p1, p2=p2, r1=r1, r2=r2,
#'                         rho=1, type= 'MatrixT', nu=3)
#' str(datlist)
#' res <- ER.RMFM(datlist$X, r_max=10,   epsELBO=1e-9, maxIter=10, verbose=FALSE, seed=1)
#' res


ER.RMFM <- function(X, r_max=10,   epsELBO=1e-9, maxIter=20, verbose=FALSE, seed=1){
  
  #library(irlba)
  message("Choose the number of factors...")
  fit <- RMFM(X, r1=r_max, r2=r_max,  epsELBO, maxIter, verbose, seed,
               cal_eigs=TRUE)
  svr_fun <- function(sv){
    n1 <- length(sv)
    sv[-n1]/ sv[-1]
  }
  svMat <- sapply(fit[c("dR", "dC")], svr_fun)
  row.names(svMat) <- 1:(r_max-1)
  rvec <- apply(svMat, 2, which.max)
  message("Estimated r1=", rvec[1], ", r2=", rvec[2]);
  return(list(rvec=rvec, svrMat=svMat))
}

