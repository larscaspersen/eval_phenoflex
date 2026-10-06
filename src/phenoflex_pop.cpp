#include <Rcpp.h>
#include <cmath> 
#include <random>
using namespace Rcpp;


namespace {
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
     if(z[i+1] >= zc && bloomindex == 0) {
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

// phenoflex_slim: compute heat accumulation using precomputed chill (x, xs, y)
List phenoflex_slim(
    NumericVector temp,
    NumericVector times,
    NumericVector x,
    NumericVector xs,
    NumericVector y,
    const double yc,
    const double zc,
    const double s1,
    const double Tb,
    const double Tu,
    const double Tc,
    const double Delta = 4,
    const int Imodel = 0,
    bool stopatzc = true,
    bool deg_celsius = true
) {
  const int N = temp.size();
  NumericVector z(N);
  z[0] = 0.0;
  
  // temperature thresholds (convert to Kelvin if needed)
  double _Tu = Tu;
  double _Tc = Tc;
  double _Tb = Tb;
  
  if (deg_celsius) {
    _Tu += 273.0;
    _Tc += 273.0;
    _Tb += 273.0;
  }
  
  int bloomindex = 0;
  
  // main loop
  for (int i = 0; i < N - 1; i++) {
    
    double ti = temp[i];
    if (deg_celsius)
      ti += 273.0;
    
    // heat increment
    if (Imodel == 0) {
      z[i + 1] = z[i] +
        P1z(ti, _Tu, _Tb, _Tc) *
        PFcn(y[i], yc, s1) *
        (times[i + 1] - times[i]);
    } else {
      z[i + 1] = z[i] +
        P2z(ti, _Tu, Delta) *
        PFcn(y[i], yc, s1) *
        (times[i + 1] - times[i]);
    }
    
    // detect bloom
    if (z[i + 1] >= zc && bloomindex == 0) {
      bloomindex = i + 2;   // Fortran-style index
      if (stopatzc)
        break;
    }
  }
  
  return List::create(
    Named("z") = z,
    Named("bloomindex") = bloomindex
  );
}

List run_forcing_experiment(NumericVector y, 
                            NumericVector z,
                              const double yc, 
                              const double zc,
                              const double s1, 
                              const int timestep,
                              const int max_days_forcing,
                              const double forcing_increment,
                              const int placeholder_fail = 999,
                              const double adjust_zc = 1){
  
  const int hours = 24; //numbers of hours forcing per day
  double y_i = y[timestep]; // get amount of chill at forcing experiment 
  const double share_effective_heat_i = PFcn(y_i, yc, s1);
  int bloom_met = placeholder_fail; //only overwrite if bloom is met

  //amount of forcing already accumulated before cutting
  double base_force = z[timestep];
  
  double eff_inc = forcing_increment * share_effective_heat_i;
  double max_force = base_force + eff_inc * (max_days_forcing * hours);
  
  //only test for budbreak when it can be reached or surprassed
  //after maximum days of forcing
  if(max_force >= (zc * adjust_zc)){
    double exp_force = base_force;
    
    //iterate through one experiment, each time step add forcing
    for (int k = 0; k < max_days_forcing*hours; k++) {
      
      exp_force += eff_inc;
      if(exp_force >= (zc * adjust_zc)){
        bloom_met = k;
        break;
      }

    } //end loop days in forcing experiments
  }
  

  
  return(List::create(Named("timestep_bloom") = bloom_met,
                      Named("share_effective_heat") = share_effective_heat_i));
}


} // namespace: keep population heat scaling separate from other model kernels

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
                   const double adjust_zc_forcing_exp = 1,
                   const int seed = 12345,
                   bool stopatzc = true,
                   bool deg_celsius = true,
                   bool basic_output = true){
  
  //unwrap the vector
  NumericVector i_cut_vec = i_cut.isNull() ? NumericVector(0) : as<NumericVector>(i_cut);
  
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
    if (Imodel == 0){
      forcing_increment = P1z(forcing_temperature, Tu, Tb, Tc);
    }
    if (Imodel == 1){
      forcing_increment =  P2z(forcing_temperature, Tu, Delta);
    }
 
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
                          stopatzc && n_experiments == 0,
                          deg_celsius,
                          basic_output && !return_force_exp);
    bloom_ind[i] = pheno_out["bloomindex"];
    
    //reset the vector tracking the forcing experiment
    NumericVector exp_bloom_i = exp_bloom_base;
    Rcpp::NumericVector share_effective_heat(n_experiments);
    
    if (!basic_output || return_force_exp){
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
                                                  placeholder_fail,
                                                  adjust_zc_forcing_exp);
        
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



//slim version of population model because chill needs to be only calculated
//once for all the population members
// [[Rcpp::export]]
List PhenoFlex_pop_slim(
    NumericVector temp,
    NumericVector times,
    const NumericVector yc,
    const NumericVector zc,
    const int max_days_forcing = 50,
    Nullable<NumericVector> i_cut = R_NilValue,
    const double forcing_temperature = 23,
    const double s1 = 0.5,
    const double E0 = 3372.8,
    const double E1 = 9900.3,
    const double A0 = 6319.5,
    const double A1 = 5.939917e13,
    const double Tf = 4,
    const double slope = 1.6,
    const double Tb = 4,
    const double Tu = 26,
    const double Tc = 36,
    const double Delta = 4,
    const int Imodel = 0,
    const int placeholder_fail = 9999,
    const double adjust_zc_forcing_exp = 1,
    bool stopatzc = true,
    bool deg_celsius = true,
    bool basic_output = true
) {
  
  // unwrap i_cut
  NumericVector i_cut_vec = i_cut.isNull() ? NumericVector(0) : as<NumericVector>(i_cut);
  
  const int N = yc.size();
  const int T = times.size();
  const int E = i_cut_vec.size();
  if (N < 1 || zc.size() != N || T < 2 || temp.size() != T) {
    stop("Supply equal non-empty requirement vectors and at least two matching temperature/time values.");
  }
  for (int i = 0; i < N; ++i) {
    if (!R_finite(yc[i]) || !R_finite(zc[i]) || yc[i] <= 0 || zc[i] <= 0)
      stop("Bud requirements must be finite and positive.");
  }
  for (int t = 0; t < T; ++t) {
    if (!R_finite(temp[t]) || !R_finite(times[t]) ||
        (t > 0 && times[t] <= times[t - 1]))
      stop("Temperatures and times must be finite and times strictly increasing.");
  }
  for (int e = 0; e < E; ++e) {
    if (!R_finite(i_cut_vec[e]) || i_cut_vec[e] != std::floor(i_cut_vec[e]) ||
        i_cut_vec[e] < 0 || i_cut_vec[e] >= T)
      stop("Cutting indices must be zero-based valid input rows.");
  }
  
  bool return_force_exp = !i_cut.isNull();
  
  // forcing increment
  double forcing_increment = return_force_exp ?
  P1z(forcing_temperature, Tu, Tb, Tc) : 0.0;
  
  // ---------------------------------------------------------
  // 1. Run full model once to compute chill (x, y, xs)
  // ---------------------------------------------------------
  List full = PhenoFlex(
    temp, 
    times,
    yc[0], 
    zc[0],
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
    false, 
    deg_celsius,
    false
  );
  
  NumericVector x = full["x"];
  NumericVector y = full["y"];
  NumericVector xs = full["xs"];
  
  // ---------------------------------------------------------
  // 2. Allocate outputs
  // ---------------------------------------------------------
  NumericVector bloom_ind(N);
  NumericMatrix z_pop(basic_output ? 0 : T, basic_output ? 0 : N);
  NumericMatrix exp_bloom_pop(E, N);
  NumericMatrix z_effec_exp(E, N);
  
  NumericVector exp_bloom_i(E);
  NumericVector share_effective_heat(E);
  
  // ---------------------------------------------------------
  // 3. Loop over population members
  // ---------------------------------------------------------
  for (int i = 0; i < N; i++) {
    
    //also calculate for i==0 because bloomindex otherwise miscalculated
    List pheno_out = phenoflex_slim(
        temp, 
        times,
        x, 
        xs,
        y,
        yc[i], 
        zc[i],
        s1,
        Tb, 
        Tu, 
        Tc,
        Delta,
        Imodel,
        // Cutting experiments need heat even after the first bloom event.
        stopatzc && E == 0,
        deg_celsius
      );
    
    bloom_ind[i] = pheno_out["bloomindex"];
    
    // store z trajectory
    NumericVector z_i = pheno_out["z"];
    if (!basic_output) {
      std::copy(z_i.begin(), z_i.end(), z_pop.begin() + i * T);
    }
    
    // ---------------------------------------------------------
    // 4. Forcing experiments
    // ---------------------------------------------------------
    if (return_force_exp) {
      
      std::fill(exp_bloom_i.begin(), exp_bloom_i.end(), placeholder_fail);
      
      for (int j = 0; j < E; j++) {
        
        List forcing_out = run_forcing_experiment(
          y, z_i,
          yc[i], zc[i],
                   s1,
                   i_cut_vec[j],
                            max_days_forcing,
                            forcing_increment,
                            placeholder_fail,
                            adjust_zc_forcing_exp
        );
        
        exp_bloom_i[j] = forcing_out["timestep_bloom"];
        share_effective_heat[j] = forcing_out["share_effective_heat"];
      }
      
      // write results
      std::copy(exp_bloom_i.begin(), exp_bloom_i.end(),
                exp_bloom_pop.begin() + i * E);
      
      std::copy(share_effective_heat.begin(), share_effective_heat.end(),
                z_effec_exp.begin() + i * E);
    }
  }
  
  // ---------------------------------------------------------
  // 5. Return output
  // ---------------------------------------------------------
  if (basic_output) {
    return List::create(Named("bloomindex") = bloom_ind);
  }
  
  return List::create(
    Named("bloomindex") = bloom_ind,
    Named("x") = x,
    Named("y") = y,
    Named("z") = z_pop,
    Named("exp") = exp_bloom_pop,
    Named("z_effec") = z_effec_exp,
    Named("exp_force_inc") = forcing_increment
  );
}
