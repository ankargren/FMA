// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List EstAllModels(arma::mat X, arma::mat Xnew, arma::mat y, arma::mat s) {

  // p: number of regressors, M: number of models, n: sample size
  // Initialize everything
  int p = s.n_cols;
  int M = s.n_rows;
  double n = X.n_rows;
  mat output(p, M);
  mat idmat = eye<mat>(p, p);
  mat etilde(n, M);
  mat eJMA(n, M);
  mat preds(Xnew.n_rows, M);
  vec SEs(M);
  vec AIC(M);
  vec BIC(M);
  vec K(M);

  // Loop over all models
  for (int ii=0; ii<M; ii++) {
    // Get the subset of variables and estimate the coefficients in the current submodel
    mat x = X.cols(find(s.row(ii) == 1));
    mat idmatsub = idmat.cols(find(s.row(ii) == 1));
    vec coefvec = inv(trans(x) * x) * trans(x) * y;
    output.col(ii) = idmatsub * coefvec;

    // Get residuals, leverage and predictions
    vec residuals = y - x * coefvec;
    etilde.col(ii) = residuals;
    vec hatvec = diagvec(x * inv(trans(x) * x) * trans(x));
    eJMA.col(ii) = etilde.col(ii) / (1 - hatvec);
    preds.col(ii) = Xnew * output.col(ii);

    // Get standard errors, number of coefficients and AIC/BIC
    double ncoef = idmatsub.n_cols;
    K[ii] = ncoef;
    SEs[ii] = sum(trans(residuals) * residuals)/(n - ncoef);
    AIC[ii] = n * log(SEs[ii]) + 2*ncoef;
    BIC[ii] = n * log(SEs[ii]) + log(n)*ncoef;
  }

  // Return a list
  return Rcpp::List::create(Rcpp::Named("SEs") = SEs,
                            Rcpp::Named("coefficients") = output, Rcpp::Named("K") = K,
                            Rcpp::Named("etilde") = etilde, Rcpp::Named("eJMA") = eJMA,
                            Rcpp::Named("AIC") = AIC, Rcpp::Named("BIC") = BIC,
                            Rcpp::Named("preds") = preds, Rcpp::Named("s") = s);
}

