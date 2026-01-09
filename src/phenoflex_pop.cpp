#include <Rcpp.h>
#include <cmath> 
#include <random>
using namespace Rcpp;

const double Pi = 3.14159265358979323846;

//GDH model
inline double P1z(const double T, const double Tu, const double Tb, const double Tc) {
  if(T >= Tb && T <= Tu) {
    return (1./2. * (1 + cos(Pi + Pi * (T - Tb)/(Tu - Tb))) );
  }
  else if(T > Tu && T <= Tc) {
    return ( (1 + cos(Pi/2. + Pi/2. * (T -Tu)/(Tc - Tu))) ); 
  }
  return (0.);
}

//alternative heat model
inline double P2z(const double T, const double Tu, const double Delta) {
  return( exp(-((T - Tu)/2./Delta)*((T - Tu)/2./Delta)) );
}

//share of conversion labile chill to heat-stable chill 
//also controls share of effective heat accumulation depending on chill
inline double PFcn(const double T, const double Tf, const double slope) {
  const double x = slope*Tf*(T-Tf)/T;
  if(x >= 17) return(1);
  else if(x <= -20) return(0);
  const double sr = exp(x);
  return( sr/(1+sr) );
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
   
   //iterate over time steps
   for(int i = 0; i < N-1; i++) {
     double ti = temp[i];
     if(deg_celsius) ti += 273.;
     
     //chill increment and destruction
     xs[i] = A0/A1 * exp(-(E0-E1)/ti);
     const double k1 = A1*exp(-E1/ti);
     
     //update pools
     x[i+1] = xs[i] - (xs[i] - x[i])*exp(-k1*(times[i+1]-times[i]));
     y[i+1] = y[i];
     
     //calculate potential heat increment
     if(Imodel == 0) {
       z[i+1] = z[i] + P1z(ti, _Tu, _Tb, _Tc) * PFcn(y[i], yc, s1)*(times[i+1]-times[i]);
     }
     else {
       z[i+1] = z[i] + P2z(ti, _Tu, Delta) * PFcn(y[i], yc, s1)*(times[i+1]-times[i]);
     }
     
     //convert labile to stable chill
     if(x[i+1] >= 1.) {
       double delta = PFcn(ti, _Tf, slope) * x[i+1];
       y[i+1] += delta;
       x[i+1] -= delta;
     }
     
     //detect if heat requirement is met
     if(z[i+1] >= zc) {
       // i+2 for Fortran index convention in R
       bloomindex = i+2;
       if(stopatzc) break;
     }
     
     //return output
   }
   if(basic_output) {
     return List::create(Named("bloomindex") = bloomindex);
   }
   return List::create(Named("x") = x, Named("y") = y, Named("z")=z, Named("xs")=xs, Named("bloomindex") = bloomindex);
 }

List run_forcing_experiment(NumericVector y, 
                            NumericVector z,
                              const double yc, 
                              const double zc,
                              const double s1, 
                              const int timestep,
                              const int max_days_forcing,
                              const double forcing_increment,
                              const int placeholder_fail = 999){
  
  const int hours = 24; //numbers of hours forcing per day
  double y_i = y[timestep]; // get amount of chill at forcing experiment 
  const double share_effective_heat_i = PFcn(y_i, yc, s1);
  int bloom_met = placeholder_fail; //only overwrite if bloom is met

  //amount of forcing already accumulated before cutting
  double base_force = z[timestep];
  
  //iterate through one experiment, each time step add forcing
  for (int k = 0; k < max_days_forcing*hours; k++) {
    //assume constant temperatures, each increment is same 
    double exp_force = base_force +  (forcing_increment* share_effective_heat_i * k);
    
    //check if zc is met
    if (exp_force >= zc){
      bloom_met = k;
      break;
    } 
  } //end loop days in forcing experiments
  
  return(List::create(Named("timestep_bloom") = bloom_met,
                      Named("share_effective_heat") = share_effective_heat_i));
}


// [[Rcpp::export]]
List PhenoFlex_pop(NumericVector temp,
                   NumericVector times,
                   const Rcpp::NumericVector yc,
                   const Rcpp::NumericVector zc,
                   const int max_days_forcing = 50,
                   Nullable<NumericVector> i_cut = R_NilValue ,
                   const int forcing_temperature = 23,
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
                   const int placeholder_fail = 9999,
                   const int seed = 12345,
                   bool stopatzc = true,
                   bool deg_celsius = true,
                   bool basic_output = true){
  
  //unwrap the vector
  NumericVector i_cut_vec(i_cut.get());
  
  const int N = yc.size();
  const int time_steps = times.size();
  const int n_experiments = i_cut_vec.size();
  double forcing_increment = 0;
  //double z_effec = 0; //share of potential heat that becomes effective (forcing experiment)
  bool return_force_exp = false;
  
  //store forcing experiment results (time steps needed to reach zc)
  NumericVector exp_bloom_base(i_cut_vec.size(), placeholder_fail);
  
  // if forcing experiment days are supplied (i_cut)
  if (i_cut.isNull() == false) { 
    return_force_exp = true; 

    //calculate the forcing potential for each day of experiment
    //assume constant conditions
    forcing_increment = P1z(forcing_temperature, Tu, Tb, Tc);
  }
  
  //where results will be saved
  Rcpp::NumericVector bloom_ind(N); //vector with time step of bloom (for each member of population)
  NumericMatrix x_pop(time_steps, N); //time step per row, population member per column
  NumericMatrix y_pop(time_steps, N); //time step per row, population member per column
  NumericMatrix z_pop(time_steps, N); //time step per row, population member per column
  NumericMatrix exp_bloom_pop(n_experiments, N); //experiment: time step after cutting zc reached

  ////////////////////
  //DEBUGGING
  NumericMatrix z_effec_exp(n_experiments, N); //share of potential heat that gets accumulated
  // END DEBUGGING
  ////////////////////
  
  //run phenoflex, return complex output
  for (int i = 0; i < N; ++i) { 
    
    //basic bloom calculation
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
    NumericVector exp_bloom_i = exp_bloom_base;
    Rcpp::NumericVector share_effective_heat(n_experiments);
    
    if (basic_output == false){
      //extract values from phenoflex list
      x_pop(_, i) = as<NumericVector>(pheno_out["x"]);
      y_pop(_, i) = as<NumericVector>(pheno_out["y"]);
      z_pop(_, i) = as<NumericVector>(pheno_out["z"]);
    }
    
    // do forcing experiment
    if(return_force_exp){
      
      //iterate over forcing experiment days
      for (int j = 0; j<n_experiments; j++) {
        
        List forcing_out = run_forcing_experiment(y_pop.column(i),
                                                  z_pop.column(i),
                                                  yc[i],
                                                  zc[i],
                                                  s1,
                                                  i_cut_vec[j],
                                                  max_days_forcing,
                                                  forcing_increment,
                                                  placeholder_fail);
        
        //save results of experiment to vectors
        share_effective_heat[j] = forcing_out["share_effective_heat"];
        exp_bloom_i[j] = forcing_out["timestep_bloom"];
        
        
        //NumericVector chill = y_pop.column(i); // extract the i-th NumericVector 
        //double val = chill[i_cut_vec[j]]; // get amount of chill at forcing experiment 
        //z_effec = PFcn(val, yc[i], s1); //calculate effectivity of forcing experiment

        ////////
        //DEBUG
        //share_effective_heat[j] = z_effec;
        ////////
        
        //NumericVector  z_vec = z_pop.column(i); //extract heat accumulated for that individual
        
        //amount of forcing already accumulated before cutting
        //double base_force = z_vec[i_cut_vec[j]];
        

        //timesteps within experiment
        //amount of forcing that was already accumulated + extra forcing from experiment * effectivity (controlled by y, yc, s1)
        //for (int k = 0; k < max_days_forcing*hours; k++) {
        //  double exp_force = base_force +  (forcing_increment* z_effec * k);
          
          //check if zc is met
        //  if (exp_force >= zc[i]){
        //    exp_bloom_i[j] = k;
        //    break;
        //  } 

        //} //end loop days in forcing experiments
        
 
      } //end loop time steps forcing experiments
    } //end experiment
    
    //save result to matrix
    exp_bloom_pop(_, i) = as<NumericVector>(exp_bloom_i);
    //std::copy(exp_bloom_i.begin(), exp_bloom_i.end(), exp_bloom_pop.begin() + i * n_experiments);
    
    ////////
    //DEBUG
    z_effec_exp(_, i) = as<NumericVector>(share_effective_heat);
    //std::copy(share_effective_heat.begin(), share_effective_heat.end(), z_effec_exp.begin() + i * n_experiments);
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
                      Named("exp_force_inc") = forcing_increment));
}
