% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gendata_rmfm.R
\name{gendata_rmfm}
\alias{gendata_rmfm}
\title{Generate simulated data}
\usage{
gendata_rmfm(
  i = 1,
  Tt = 100,
  p1 = 50,
  p2 = 40,
  r1 = 4,
  r2 = 3,
  rho = 0.01,
  type = c("MatrixT", "MatrixN"),
  nu = 1
)
}
\arguments{
\item{i}{a positive integer, the random seed for reproducibility of data generation process.}

\item{Tt}{a positive integer, specify the sample size.}

\item{p1}{a positive integer, specify the row dimension of the observed matrix.}

\item{p2}{a positive integer, specify the column dimension of the observed matrix.}

\item{r1}{a positive integer, specify the number of row factors; default as 4}

\item{r2}{a positive integer, specify the number of column factors; default as 3.}

\item{rho}{a positive real, specify the signal strength of factor matrices.}

\item{type}{a string, specify the type of error matrix, default as \code{type='MatrixN'}; supportint matrix t distribution 'MatrixT' and matrix normal distribution 'MatrixN'.}

\item{nu}{a positive integer, specify the degree freedom of the matrix t distribution when \code{type='MatrixT'}.}
}
\value{
return a list including the following components:
\itemize{
  \item \code{X} - p1* p2*T array, which is the observed  matrix from each individual, where T is the sample size.
  \item \code{CC} - p1* p2*T array, which is the common component matrix for each individual.
  \item \code{F0} - r1* r2*T array, which is the generated factor  matrix for each individual, where T is the sample size.
  \item \code{R0} - a p1-by-r1 matrix, the row loading matrix.
  \item \code{C0} - a p2-by-r2 matrix, the column loading matrix.
  \item \code{mu0} - a p1-by-p2 matrix, the mean matrix.
}
}
\description{
Generate simulated data from robust matrix factor models
}
\examples{
r1 <- 4; r2 <- 3;
Tt <- 100; type <- 'MatrixT'
p1 <- 100; p2 <- 50
datlist <- gendata_rmfm(i = 1,  Tt = Tt,p1 =p1, p2=p2, r1=r1, r2=r2,
                         rho=0.01, type=type, nu=1)
str(datlist)
}
