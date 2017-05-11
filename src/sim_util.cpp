#include "sim.h"
//
/***********************************************************************************/
/* myacf: autocovariance function                                                  */
/***********************************************************************************/
//
arma::vec myacf(arma::vec const &ts,
                int lag){
  double    mean = arma::mean(ts);  int N = ts.n_elem;
  arma::vec out(lag+1);
  for(int l=0; l<=lag; l++){
    double autocv=0.0;
    for(int i=0; i<(N-l); i++){
      autocv+=((ts[i]-mean)*(ts[i+l]-mean));
    }
    out[l] = (1.0/(N)) * autocv;
  }
  return out;
}
//
/***********************************************************************************/
/* specParzen: Parzen window smoothed PSD                                          */
/***********************************************************************************/
//[[Rcpp::export]]
arma::vec specParzen(arma::vec const &ts,
                     int lag,
                     int const &maxf,
                     int const &outn){
  // maximum lag as N-1;
  lag = (int) (lag<ts.n_elem)? lag : (ts.n_elem-1);
  arma::vec xcovs = myacf(ts, lag);                     //acf
  arma::vec win   = linspace(0,1,(lag+1));              //parzen win
  // Parzen window of lag+1;
  for(int i=0; i<=lag; i++){
    win[i] = (double) (win[i] < .5)? (1.0 - 6.0 *
        std::pow(win[i],2.0) * (1.0 - win[i])):
        (2.0*std::pow((1-win[i]), 3.0));
  }
  xcovs %= win;                                         //smoothed acf
  double up = ((double) maxf)/ts.n_elem;
  arma::vec spec = arma::linspace((up/outn),up,outn);
  for(int i=0; i<outn; i++){
    double sum=xcovs[0];
    for(int ii=1; ii<=lag; ii++){
      sum += 2.0 * cos (2.0*spec[i]*datum::pi*ii) * xcovs[ii];
    }
    spec[i] = sum;
  }
return spec;
}
//
/***********************************************************************************/
/* MA: Moving Average                                                              */
/***********************************************************************************/
//
arma::cube MA(arma::cube spec,
              int const &win,
              int const &o){
  if(win == 1) return spec;
  int end = win-1;            int start = 0;
  int ns = spec.n_slices;     int ne = 0;
  arma::cube out(size(spec)); out.fill(0.0);
  while(end<ns){
    for(int ii=0; ii<win; ii++){
      out.slice(ne) += spec.slice(start+ii) / win;
    }
    start += (int) (win - o);
    end   += (int) (win - o);
    ne++;
  }
  out.shed_slices(ne, ns-1);
  return out;
}
//
/***********************************************************************************/
/* SpecSim: similarity between spectral densities                                  */
/***********************************************************************************/
//
//[[Rcpp::export]]
arma::cube SpecSim(arma::cube const &ts,
              int lag, int const &wn,
              int const &win, int const &overlap,
              int const &specN){
  // Estiamte spectral array
  int nc = ts.n_cols;     int ns = ts.n_slices;
  arma::cube segspec(specN, nc, ns);
  for(int seg=0; seg<ns; seg++){
    for(int chan=0; chan<nc; chan++){
      arma::vec vts = arma::conv_to<vec>::from(ts.slice(seg).col(chan));
      arma::vec vsp = specParzen(vts, lag, wn, specN);
      vsp /= arma::accu(vsp);
      segspec.slice(seg).col(chan) = vsp;
    }
  }
  // Moving average
  arma::cube epspec = MA(segspec, win, overlap);
  int nep = epspec.n_slices;
  // distance
  arma::cube epdis(nc,nc,nep);
  for(int ne=0; ne<nep; ne++){
    for(int c1=0; c1<nc; c1++){
      for(int c2=c1; c2<nc; c2++){
        epdis(c1,c2,ne) = (c1==c2)? (double) 0.0:
        (arma::accu(arma::min(epspec.slice(ne).col(c1),
                              epspec.slice(ne).col(c2))));
        epdis(c2,c1,ne) = epdis(c1,c2,ne);
      }
    }
  }
  return epdis;
}
//
/***********************************************************************************/
/* EigLap: Eigen Laplacian transformation                                          */
/***********************************************************************************/
//
//[[Rcpp::export]]
arma::cube EigLap(
    arma::cube  const & data,
    int         const &D,
    bool        const &normal){
  int nc = data.n_rows;     int ne = data.n_slices;
  arma::cube eigdata(D, nc, ne);
  for(int ie=0; ie<ne; ie++){
    arma::vec eigval; arma::mat eigvec;
    if(normal){ // with normalized Laplacian;
      arma::colvec rsum = sum(data.slice(ie), 1);
      arma::mat    nmat = arma::diagmat(arma::conv_to<vec>::from(pow(rsum, -0.5)));
      arma::mat    ndat = nmat * data.slice(ie) * nmat;
      eig_sym(eigval, eigvec, ndat);
      if(D < nc) eigvec.shed_cols(0,((int) (nc-D-1)));
      for(int ir=0;ir<nc;ir++){
        double cc = std::pow(arma::accu(pow(eigvec.row(ir),2.0)), -0.5);
        eigvec.row(ir) *= cc;
      }
    } else {    //Unnormalized Lap
      eig_sym(eigval, eigvec, data.slice(ie));
      if(D < nc) eigvec.shed_cols(0,((int) (nc-D-1)));
    }
    eigdata.slice(ie) = eigvec.t();
  }
  return eigdata;
}
//
/***********************************************************************************/
/* EigLapSph: Eigen Laplacian on Spherical coordinates                             */
/***********************************************************************************/
//
//[[Rcpp::export]]
arma::cube EigLapSph(
    arma::cube  const & data,
    int         const &D
  ){
  int nc = data.n_rows;     int ne = data.n_slices;
  arma::cube eigdata(D, nc, ne);
  for(int ie=0; ie<ne; ie++){                           //each slice
    arma::vec eigval; arma::mat eigvec;
    arma::colvec rsum = sum(data.slice(ie), 1);
    arma::mat    nmat = arma::diagmat(arma::conv_to<vec>::from(pow(rsum, -0.5)));
    arma::mat    ndat = nmat * data.slice(ie) * nmat;   //normalized datamat
    eig_sym(eigval, eigvec, ndat);
    if(D < nc) eigvec.shed_cols(0,((int) (nc-D-1)));    //preserving D-dim eigvecs
    for(int ir=0;ir<nc;ir++){                           //each channel(row) ir
      //  radius
      eigdata(0,ir,ie) = std::pow(arma::accu(pow(eigvec.row(ir),2.0)), -0.5);
      //  angles
      for(int ic=1;ic<D; ic++){
        double cc = std::pow(arma::accu(arma::pow(eigvec(ir, span(ic-1, D-1)), 2.0)),-0.5);
        eigdata(ic,ir,ie) = std::acos ((double) (eigvec(ir,(ic-1))/cc));
      }
      eigdata((D-1),ir,ie) = (eigvec(ir,(D-1))>=0)? eigdata((D-1),ir,ie):
        2*datum::pi - eigdata((D-1),ir,ie);
    }
    eigdata.slice(ie) = eigvec.t();
  }
  return eigdata;
}
