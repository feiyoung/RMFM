#' Generate simulated data
#' @description Generate simulated data from robust matrix factor models
#' @param i a positive integer, the random seed for reproducibility of data generation process.
#' @param Tt a positive integer, specify the sample size.
#' @param p1 a positive integer, specify the row dimension of the observed matrix.
#' @param p2 a positive integer, specify the column dimension of the observed matrix.
#' @param r1 a positive integer, specify the number of row factors; default as 4
#' @param r2 a positive integer, specify the number of column factors; default as 3.
#' @param rho a positive real, specify the signal strength of factor matrices.
#' @param type a string, specify the type of error matrix, default as \code{type='MatrixN'}; supportint matrix t distribution 'MatrixT' and matrix normal distribution 'MatrixN'.
#' @param nu a positive integer, specify the degree freedom of the matrix t distribution when \code{type='MatrixT'}.
#' @return return a list including the following components:
#' \itemize{
#'   \item \code{X} - p1* p2*T array, which is the observed  matrix from each individual, where T is the sample size.
#'   \item \code{CC} - p1* p2*T array, which is the common component matrix for each individual.
#'   \item \code{F0} - r1* r2*T array, which is the generated factor  matrix for each individual, where T is the sample size.
#'   \item \code{R0} - a p1-by-r1 matrix, the row loading matrix.
#'   \item \code{C0} - a p2-by-r2 matrix, the column loading matrix.
#'   \item \code{mu0} - a p1-by-p2 matrix, the mean matrix.
#' }
#' @importFrom LaplacesDemon rmatrixnorm
#' @importFrom MixMatrix  rmatrixt
#' @export
#'
#' @examples
#' r1 <- 4; r2 <- 3;
#' Tt <- 100; type <- 'MatrixT'
#' p1 <- 100; p2 <- 50
#' datlist <- gendata_rmfm(i = 1,  Tt = Tt,p1 =p1, p2=p2, r1=r1, r2=r2,
#'                          rho=0.01, type=type, nu=1)
#' str(datlist)


gendata_rmfm <- function(i=1, Tt=100,  p1 =50, p2=40, r1=4, r2=3,
                            rho=0.01, type=c("MatrixT", "MatrixN"), nu=1){
 
  type <- match.arg(type)
  # if(type=='MatrixN'){
  #   library(LaplacesDemon)
  # }else if(type=='MatrixT'){
  #   library(MixMatrix)
  # }
  
  
  # install.packages('MixMatrix')
  set.seed(i) # 50
  Ff <- array(dim=c(r1,r2, Tt))
  ## generate F, Y and X;
  mu <- matrix(rnorm(p1*p2, mean = 1), p1, p2)
  R0 <- mat2orth(matrix(rnorm(p1*r1), p1, r1))
  C0 <- mat2orth(matrix(rnorm(p2*r2), p2, r2))
  #Lambda10 <- rep(1, p1)/p1; Lambda20 <- rep(0.5, p2) /p2
  Y <- array(dim=c(p1,p2, Tt))
  CP <- Y
  
  Eps <- array(0,dim=c(r1,r2, Tt))
  Et <- array(0,dim=c(p1,p2, Tt))
  Ut <- array(0, dim=c(p1,p2, Tt))
  phi <- 0.6; psi <- 0.4
  if(type=='MatrixN'){
    Eps[,,1] <- LaplacesDemon::rmatrixnorm(matrix(0, r1, r2), diag(rep(1,r1)), diag(rep(1, r2)))
    Ff[,,1] <- sqrt(1-phi^2) * LaplacesDemon::rmatrixnorm(matrix(0, r1, r2), diag(rep(1,r1)), diag(rep(1, r2)))
    Ut[,,1] <- LaplacesDemon::rmatrixnorm(matrix(0, p1, p2), diag(rep(0.5,p1)), diag(rep(0.3, p2)))
    Et[,,1] <- sqrt(1-psi^2) * LaplacesDemon::rmatrixnorm(matrix(0, p1, p2), diag(rep(0.5,p1)), diag(rep(0.3, p2)))
  }else if(type == 'MatrixT'){
    Eps[,,1] <- MixMatrix::rmatrixt(n=1, mean=matrix(0, r1, r2), L=diag(rep(1,r1)), R=diag(rep(1, r2)), df=nu)
    Ff[,,1] <- sqrt(1-phi^2) *  MixMatrix::rmatrixt(n=1, mean=matrix(0, r1, r2), L=diag(rep(1,r1)), R=diag(rep(1, r2)), df=nu)
    Ut[,,1] <-  MixMatrix::rmatrixt(n=1, mean=matrix(0, p1, p2), L=diag(rep(0.5,p1)), R=diag(rep(0.3, p2)), df=nu)
    Et[,,1] <- sqrt(1-psi^2) * MixMatrix::rmatrixt(n=1, mean=matrix(0, p1, p2), L=diag(rep(0.5,p1)), R=diag(rep(0.3, p2)), df=nu)
  }
  t <- 1
  CP[,,t] <- mu + rho*R0 %*% Ff[,,t] %*%  t(C0)
  Y[,,t] <-  CP[,,t] +  Et[,,t]
  
  for(t in 2:Tt){
    # t <- 90
    if(type=='MatrixN'){
      Ut[,,t] <- LaplacesDemon::rmatrixnorm(matrix(0, p1, p2), diag(rep(0.5,p1)), diag(rep(0.3, p2)))
      Eps[,,t] <- LaplacesDemon::rmatrixnorm(matrix(0, r1, r2), diag(rep(1,r1)), diag(rep(1, r2)))
    }else if(type == 'MatrixT'){
      Ut[,,t] <-  MixMatrix::rmatrixt(n=1, mean=matrix(0, p1, p2), L=diag(rep(0.5,p1)), R=diag(rep(0.3, p2)), df=nu)
      Eps[,,t] <-  MixMatrix::rmatrixt(n=1, mean=matrix(0, r1, r2), L=diag(rep(1,r1)), R=diag(rep(1, r2)), df=nu)
      #print(Eps[,,t][1,])
    }
    Et[,,t] <- psi* Et[,,t-1]  + sqrt(1-psi^2) * Ut[,,t]
    Ff[,,t] <- phi*Ff[,,t-1] + sqrt(1-phi^2) * Eps[,,t]
    CP[,,t] <- mu +  rho*R0 %*%Ff[,,t] %*%  t(C0)
    Y[,,t] <-  CP[,,t] +  Et[,,t]
  }
  return(list(X=Y, CC=CP, F0=Ff, R0=R0, C0=C0, mu0=mu))
  
}
