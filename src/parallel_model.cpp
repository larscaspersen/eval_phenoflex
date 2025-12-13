#include <Rcpp.h>
#include <cmath>        
using namespace Rcpp;

const double Pi = 3.14159265358979323846;

inline double P1z(const double T, const double Tu, const double Tb, const double Tc) {
  if(T >= Tb && T <= Tu) {
    return ((Tu - Tb)/2. * (1 + cos(Pi + Pi * (T - Tb)/(Tu - Tb))) );
  }
  else if(T > Tu && T <= Tc) {
    return ( (Tu - Tb) * (1 + cos(Pi/2. + Pi/2. * (T -Tu)/(Tc - Tu))) ); 
  }
  return (0.);
}

inline double Com(const double Kmin, const double yc, const double y) {
  if(y >= yc){
    return(1.);
  }
  
  return(Kmin + ((1-Kmin)/yc) * y);
  
}


inline double PFcn(const double T, const double Tf, const double slope) {
  const double x = slope*Tf*(T-Tf)/T;
  if(x >= 17) return(1);
  else if(x <= -20) return(0);
  const double sr = exp(x);
  return( sr/(1+sr) );
}

//' @title parallel_model
//' @description Parallel model, combining dynamic model for chill accumulation and the GDH model
//'
//' @param yc numeric. Critical value defining end of chill accumulation
//' @param zc numeric. Critical value of z determining the end of heat accumulation
//' @param kmin numeric. Share of buds that can flower without receiving chill. 
//' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
//' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
//' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
//' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
//' @param slope numeric. Slope parameter for sigmoidal function
//' @param Tf numeric. Transition temperature (in degree Kelvin) for the sigmoidal function
//' @param Tb numeric. GDH base temperature (lower threshold) 
//' @param Tu numeric. GDH optimal temperature 
//' @param Tc numeric. GDH upper temperature (upper threshold)
//' @param Delta numeric. Width of Gaussian heat accumulation model
//' @param stopatzc boolean. If `TRUE`, the PhenoFlex is applied until the end of the temperature series. Default is to stop once the value zc has been reached.
//' @param deg_celsius boolean. If set `TRUE` function assumes degree celsius temperature parameters, otherwise kelvin.
//' @param basic_output boolean. If `TRUE`, only the bloomindex is returned as a named element of the return list.
//' @useDynLib evalpheno
//' @author Lars Caspersen <lcaspers@uni-bonn.de>
//' @return
//' A list is returned with named element `bloomindex`, which is the index at which blooming occurs. When `basic_output=FALSE` also `x`, `y`, `z` and `xs` are
//' returned as named element of this list, which are numeric vectors of the same length as the input vector `temp` containing the hourly temperatures.
//' @examples
//' data(KA_weather)
//' hourtemps <- stack_hourly_temps(KA_weather, latitude=50.4)
//' iSeason <- genSeason(hourtemps, years=c(2009))
//' zc <- 6000
//' yc <- 40
//' kmin <- 0.1
//' x <- parallel_model(temp=hourtemps$hourtemps$Temp[iSeason[[1]]],
//'                times=c(1: length(hourtemps$hourtemps$Temp[iSeason[[1]]])),
//'                zc=zc, stopatzc=TRUE, yc=yc, kmin = kmin, basic_output=FALSE)
//' DBreakDay <- x$bloomindex
//' ii <- c(1:DBreakDay)
//' plot(x=ii, y=x$z[ii], xlab="Hour Index", ylab="z", col="red", type="l")
//' abline(h=zc, lty=2)
//' plot(x=ii, y=x$y[ii], xlab="Hour Index", ylab="y", col="red", type="l")
//' abline(h=yc, lty=2)
//' @export
// [[Rcpp::export]]
 List parallel_model(NumericVector temp,
               NumericVector times,
               const double yc=40,
               const double zc=1119,
               const double kmin=0.1,
               const double A0=6319.5,
               const double A1=5.939917e13,
               const double E0=3372.8,
               const double E1=9900.3,
               const double slope=1.6,
               const double Tf=4,
               const double Tu=25,
               const double Tb=4,
               const double Tc=36,
               const double Delta=4,
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
    if(y[i] < yc){
     x[i+1] = xs[i] - (xs[i] - x[i])*exp(-k1*(times[i+1]-times[i]));
     y[i+1] = y[i];
    } 
    else {
     x[i+1] = x[i];
     y[i+1] = y[i];
    }
    
    if(x[i+1] >= 1.) {
     double delta = PFcn(ti, _Tf, slope) * x[i+1];
     y[i+1] += delta;
     x[i+1] -= delta;
    }
    
   // calculate heat increment and add to z
   z[i+1] = z[i] + (Com(kmin, yc, y[i+1]) * P1z(ti, _Tu, _Tb, _Tc));
   
 
 // if z_crit calculated (and not zero as default)
 // check if it is reached, if so return bloom index
 if(z[i+1] >= zc){
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
