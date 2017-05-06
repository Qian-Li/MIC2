#include "sim.h"
#include "MIC.h"

// arma::mat TVD()
// {
//   arma::vec X  = linspace<vec> (0, datum::pi/2.0, 1000);
//   arma::vec Y1 = sin(X);  arma::vec Y2 = cos(X);
//   arma::vec D = arma::min(Y1, Y2);
//
//   return trapz(X,D);
// }

// // [[Rcpp::export]]
// void gmm_test(){
//   uword d = 3;       // dimensionality
//   uword N = 12;   // number of vectors
//
//   mat data(d, N, fill::zeros);
//
//   vec mean0 = linspace<vec>(1,d,d);
//   vec mean1 = mean0 + 2;
//
//   uword i = 0;
//
//   while(i < N)
//   {
//     if(i < N)  { data.col(i) = mean0 + randn<vec>(d); ++i; }
//     if(i < N)  { data.col(i) = mean0 + randn<vec>(d); ++i; }
//     if(i < N)  { data.col(i) = mean1 + randn<vec>(d); ++i; }
//   }
//
//
//   // model the data as a GMM with 2 Gaussians
//
//   gmm_diag model;
//
//   bool status = model.learn(data, 2, maha_dist, random_subset, 10, 10, 1e-10, true);
//
//   if(status == false)
//   {
//     cout << "learning failed" << endl;
//   }
//   rowvec wts(2, fill::zeros); wts[1] = 1;
//   model.means.print("means:");
//   model.dcovs.print("stds:");
//   model.set_hefts(wts);
//   wts.print("wt1");
//   // double  scalar_likelihood = model.log_p( data.col(0)    );
//   rowvec     set_likelihood1 = model.log_p( data.cols(0,9) ,1);
//   rowvec     overall1 = model.log_p(data.cols(0,9));
//   set_likelihood1.print("log-like1"); overall1.print("overall prob:");
//   wts.fill(1.0/2);
//   wts.print("wt2:");
//   model.set_hefts(wts);
//   rowvec     set_likelihood2 = model.log_p( data.cols(0,9) ,0);
//   rowvec     overall2 = model.log_p(data.cols(0,9));
//   set_likelihood2.print("log-like2"); overall2.print("overall prob:");
//   // return(data);
//
//   // // Testing on mrmultinom;
//   // arma::mat llm = randu<mat>(4,10) - 1000;
//   // arma::mat sampe=zeros(size(llm));
//   // for(int i=0; i<100000; i++) sampe += mrmultinom(llm);
//   // sampe /= 100000.0;
//   // for(int nc = 0; nc<10; nc++) llm.col(nc) = arma::exp(llm.col(nc)+1000)/accu(arma::exp(llm.col(nc)+1000));
//   // llm.print("like-mat");
//   // sampe.print("sample");
// }
//
// // Testing1:
// //  data <- gmm_test()
// //  log(dmvnorm(t(data), mean = c(3.0020,4.4426,5.3218), sigma = diag(c(0.6088,0.1630,0.1845))))
// // Testing2:
// //  Sys.time()->start; gmm_test(); Sys.time()-start
