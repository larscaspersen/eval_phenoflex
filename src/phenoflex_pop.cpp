#include <Rcpp.h>
#include <cmath> 
#include <random>
using namespace Rcpp;

const double Pi = 3.14159265358979323846;

inline double P1z(const double T, const double Tu, const double Tb, const double Tc) {
  if(T >= Tb && T <= Tu) {
    return (1./2. * (1 + cos(Pi + Pi * (T - Tb)/(Tu - Tb))) );
  }
  else if(T > Tu && T <= Tc) {
    return ( (1 + cos(Pi/2. + Pi/2. * (T -Tu)/(Tc - Tu))) ); 
  }
  return (0.);
}

inline double P2z(const double T, const double Tu, const double Delta) {
  return( exp(-((T - Tu)/2./Delta)*((T - Tu)/2./Delta)) );
}

// [[Rcpp::export]]
inline double PFcn(const double T, const double Tf, const double slope) {
  const double x = slope*Tf*(T-Tf)/T;
  if(x >= 17) return(1);
  else if(x <= -20) return(0);
  const double sr = exp(x);
  return( sr/(1+sr) );
}

// draw values from normal distribution
std::vector<double> gen_pop(const double mean, const double sd, const int n, const int seed = 12345) {   
  std::mt19937 gen(seed); // generator

    // Normal distribution with given mean and sd 
  std::normal_distribution<> dist(mean, sd);
  
  // Vector to hold results 
  std::vector<double> values; 
  values.reserve(n); // reserve memory for efficiency
  
  // Draw n samples 
  for (int i = 0; i < n; ++i) { 
    values.push_back(dist(gen));
    }
  return values;
}

List PhenoFlex(NumericVector temp,
                NumericVector times,
                const double yc=40,
                const double zc=190,
                const double s1=0.5,
                const double E0=3372.8,
                const double E1=9900.3,
                const double A0=6319.5,
                const double A1=5.939917e13,
                const double Tf=4,
                const double slope=1.6,
                const double Tb=4,
                const double Tu=26,
                const double Tc=36,
                const double Delta=4,
                const int Imodel=0,
                bool stopatzc = true,
                bool deg_celsius = true,
                bool basic_output = true) {
   
   const int N = temp.size();
   Rcpp::NumericVector x(N);
   Rcpp::NumericVector y(N);
   Rcpp::NumericVector z(N);
   Rcpp::NumericVector xs(N);
   x[0] = 0.;
   y[0] = 0.;
   z[0] = 0.;
   double _Tf = Tf;
   double _Tu = Tu;
   double _Tc = Tc;
   double _Tb = Tb;
   if(deg_celsius) {
     _Tf += 273.;
     _Tu += 273.;
     _Tc += 273.;
     _Tb += 273.;
   }
   int bloomindex = 0;
   for(int i = 0; i < N-1; i++) {
     double ti = temp[i];
     if(deg_celsius) ti += 273.;
     xs[i] = A0/A1 * exp(-(E0-E1)/ti);
     const double k1 = A1*exp(-E1/ti);
     x[i+1] = xs[i] - (xs[i] - x[i])*exp(-k1*(times[i+1]-times[i]));
     y[i+1] = y[i];
     if(Imodel == 0) {
       z[i+1] = z[i] + P1z(ti, _Tu, _Tb, _Tc) * PFcn(y[i], yc, s1)*(times[i+1]-times[i]);
     }
     else {
       z[i+1] = z[i] + P2z(ti, _Tu, Delta) * PFcn(y[i], yc, s1)*(times[i+1]-times[i]);
     }
     if(x[i+1] >= 1.) {
       double delta = PFcn(ti, _Tf, slope) * x[i+1];
       y[i+1] += delta;
       x[i+1] -= delta;
     }
     if(z[i+1] >= zc) {
       // i+2 for Fortran index convention in R
       bloomindex = i+2;
       if(stopatzc) break;
     }
   }
   if(basic_output) {
     return List::create(Named("bloomindex") = bloomindex);
   }
   return List::create(Named("x") = x, Named("y") = y, Named("z")=z, Named("xs")=xs, Named("bloomindex") = bloomindex);
 }
// [[Rcpp::export]]
List PhenoFlex_pop(NumericVector temp,
                   NumericVector times,
                   const Rcpp::NumericVector yc,
                   const Rcpp::NumericVector zc,
                   const int max_days_forcing = 50,
                   Nullable<NumericVector> i_cut = R_NilValue ,
                   const int force_temp = 23,
                   const double s1=0.5,
                   const double E0=3372.8,
                   const double E1=9900.3,
                   const double A0=6319.5,
                   const double A1=5.939917e13,
                   const double Tf=4,
                   const double slope=1.6,
                   const double Tb=4,
                   const double Tu=26,
                   const double Tc=36,
                   const double Delta=4,
                   const int Imodel=0,
                   const int seed = 12345,
                   bool stopatzc = true,
                   bool deg_celsius = true,
                   bool basic_output = true){
  
  // draw population of yc and zc
  //std::vector<double> yc_pop = gen_pop(yc, yc_sd, n, seed);
  //std::vector<double> zc_pop = gen_pop(zc, zc_sd, n, seed+1);
  
  bool return_force_exp = false;
  //unwrap the vector
  NumericVector i_cut_vec(i_cut.get());
  
  //vector that keeps track how many time steps after cutting bloom was
  //reached in experiment
  //will be saved later to the exp_bloom_pop list
  //-999 is placeholder as failure to reach zc
  const std::vector<int> exp_bloom_base(i_cut_vec.size(), -999);
  
  
  double force_inc = 0;
  
  double z_effec = 0; //variable for forcing experiment, how much heat increment actually accumulated
  
  // override basic_output depending on i_cut 
  if (i_cut.isNull() == false) { 
    return_force_exp = true; 

    //calculate the forcing potential for each day of experiment
    //assume constant conditions
    force_inc = P1z(force_temp, Tu, Tb, Tc);
    
  }
  
  //in case both sds are 1 and we effectively only sample one value
  //run the calculations only once
  int n_run = yc.size();

  Rcpp::NumericVector bloom_ind(n_run);
  NumericMatrix x_pop(times.size(), n_run);
  NumericMatrix y_pop(times.size(), n_run);
  NumericMatrix z_pop(times.size(), n_run);
  NumericMatrix exp_bloom_pop(i_cut_vec.size(),n_run); //list with the forcing of the individual cutting days
  //vector with as many positions as experiments
  // entry says how many time steps after cutting the zc was reached
  //-999 in case it was not met
  
  
  ////////////////////
  //DEBUGGING
  
  NumericMatrix z_effec_exp(i_cut_vec.size(), n_run);

  // END DEBUGGING
  ////////////////////
  

  //run phenoflex, return complex output
  for (int i = 0; i < n_run; ++i) { 
    
    List pheno_out = PhenoFlex(temp,
                          times,
                          yc[i],
                          zc[i],
                          s1,
                          E0,
                          E1,
                          A0,
                          A1,
                          Tf,
                          slope,
                          Tb,
                          Tu,
                          Tc,
                          Delta,
                          Imodel,
                          stopatzc,
                          deg_celsius,
                          basic_output);
    bloom_ind[i] = pheno_out["bloomindex"];
    
    //reset the vector tracking the forcing experiment
    std::vector<int> exp_bloom_i = exp_bloom_base;
    Rcpp::NumericVector z_effec_vec(i_cut_vec.size());
    
    if (basic_output == false){
      //extract values from phenoflex list
      NumericVector x_int = pheno_out["x"];
      NumericVector y_int = pheno_out["y"];
      NumericVector z_int = pheno_out["z"];

      //save results to matrix
      std::copy(x_int.begin(), x_int.end(), x_pop.begin() + i * times.size());
      std::copy(y_int.begin(), y_int.end(), y_pop.begin() + i * times.size());
      std::copy(z_int.begin(), z_int.end(), z_pop.begin() + i * times.size());
      

    }
    
    
    // do forcing experiment
    if(return_force_exp){
      
      //iterate over forcing experiment days
      for (int j = 0; j<i_cut_vec.size(); j++) {
        
        NumericVector chill = y_pop.column(i); // extract the i-th NumericVector 
        double val = chill[i_cut_vec[j]]; // get amount of chill at forcing experiment 
        z_effec = PFcn(val, yc[i], s1); //calculate effectivity of forcing experiment

        ////////
        //DEBUG
        z_effec_vec[j] = z_effec;
        ////////
        
        NumericVector  z_vec = z_pop.column(i); //extract heat accumulated for that individual
        
        //amount of forcing already accumulated before cutting
        double base_force = z_vec[i_cut_vec[j]];
        

        //timesteps within experiment
        //amount of forcing that was already accumulated + extra forcing from experiment * effectivity (controlled by y, yc, s1)
        for (int k = 0; k < max_days_forcing*24; k++) {
          double exp_force = base_force +  (force_inc* z_effec * k);
          
          //check if zc is met
          if (exp_force >= zc[i]){
            exp_bloom_i[j] = k;
            break;
          } 

        } //end loop days in forcing experiments
        
 
      } //end loop time steps forcing experiments
    } //end experiment
    
    //save result to matrix
    std::copy(exp_bloom_i.begin(), exp_bloom_i.end(), exp_bloom_pop.begin() + i * i_cut_vec.size());
    
    ////////
    //DEBUG
    std::copy(z_effec_vec.begin(), z_effec_vec.end(), z_effec_exp.begin() + i * i_cut_vec.size());
    //z_effec_list[i] = z_effec_vec;
    ////////
    
  } //end loop population members
  
  if (basic_output) {
    return(List::create(Named("bloomindex") = bloom_ind));
  }
  return(List::create(Named("bloomindex") = bloom_ind,
                      Named("x") = x_pop,
                      Named("y") = y_pop,
                      Named("z") = z_pop,
                      Named("exp") = exp_bloom_pop,
                      Named("z_effec") = z_effec_exp,
                      Named("exp_force_inc") = force_inc));
}
