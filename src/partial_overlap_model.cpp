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

inline double P2z(const double T, const double Tu, const double Delta) {
  return( exp(-((T - Tu)/2./Delta)*((T - Tu)/2./Delta)) );
}

inline double PFcn(const double T, const double Tf, const double slope) {
  const double x = slope*Tf*(T-Tf)/T;
  if(x >= 17) return(1);
  else if(x <= -20) return(0);
  const double sr = exp(x);
  return( sr/(1+sr) );
}

//' @title partial_overlap_model
 //' @description Partial Overlap Model, combining the dynamic model for chill accumulation and the GDH model
 //'
 //' @param yc numeric, Critical value defining end of chill accumulation
 //' @param b1 numeric. Heat requirement at minimum chilling
 //' @param b2 numeric. Heat requiement at maximum chilling
 //' @param b3 numeric. Scales the compensation from heat and chill requirement. Low b3 leads to linear compensation, high b3 to much heat needed to compensate chill.
 //' @param ol numeric. Controls how much longer chill accumulates after reaching yc. Chill accumulation stops when share of heat requirement is surprassed.
 //' @param A0 numeric. Parameter \eqn{A_0}{A0} of the dynamic model
 //' @param A1 numeric. Parameter \eqn{A_1}{A1} of the dynamic model
 //' @param E0 numeric. Parameter \eqn{E_0}{E0} of the dynamic model
 //' @param E1 numeric. Parameter \eqn{E_1}{E1} of the dynamic model
 //' @param slope numeric. Slope parameter for sigmoidal function
 //' @param Tf numeric. Transition temperature (in degree Kelvin) for the sigmoidal functionr
 //' @param Tb numeric. GDH base temperature (lower threshold)  //' @param Tc numeric. GDH upper temperature (upper threshold)
 //' @param Tu numeric. GDH optimal temperature 
 //' @param Tc numeric. GDH critical temperature 
 //' @param Delta numeric. Width of Gaussian heat accumulation model
 //' @param stopatzc boolean. If `TRUE`, the PhenoFlex is applied until the end of the temperature series. Default is to stop once the value zc has been reached.
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
 //' yc <- 40
 //' b1 <- 1119
 //' b2 <- 8677
 //' b3 <- 0.01119
 //' ol <- 0.75
 //' x <- po_model(temp=hourtemps$hourtemps$Temp[iSeason[[1]]],
 //'                times=c(1: length(hourtemps$hourtemps$Temp[iSeason[[1]]])),
 //'                yc=yc, b1 = b1, b2=b2, b3=b3, ol = ol, stopatzc=TRUE, basic_output=FALSE)
 //' DBreakDay <- x$bloomindex
 //' ii <- c(1:DBreakDay)
 //' plot(x=ii, y=x$z[ii], xlab="Hour Index", ylab="z", col="red", type="l")
 //' abline(h=zc, lty=2)
 //' plot(x=ii, y=x$y[ii], xlab="Hour Index", ylab="y", col="red", type="l")
 //' abline(h=yc, lty=2)
 //' @export
 // [[Rcpp::export]]
 List po_model(NumericVector temp,
                NumericVector times,
                const double yc=40,
                const double b1=1119,
                const double b2=8677,
                const double b3=0.01119,
                const double ol=0.75,
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
   const double z_chill_stop = b1 * ol;
   double z_crit = 0;
   int bloomindex = 0;
   for(int i = 0; i < N-1; i++) {
     double ti = temp[i];
     if(deg_celsius) ti += 273.;
     xs[i] = A0/A1 * exp(-(E0-E1)/ti);
     
     //' only accumulate chill when the overlap of chill and heat is not exhausted
     //' otherwise stop accumulation
     const double k1 = A1*exp(-E1/ti);
     if(z[i] < z_chill_stop){
       x[i+1] = xs[i] - (xs[i] - x[i])*exp(-k1*(times[i+1]-times[i]));
       y[i+1] = y[i];
     } 
     else {
       x[i+1] = x[i];
       y[i+1] = y[i];
     }
     
     //' in case the increment is larger than one, trunk the increment by max delta
     if(x[i+1] >= 1.) {
       double delta = PFcn(ti, _Tf, slope) * x[i+1];
       y[i+1] += delta;
       x[i+1] -= delta;
     }

     
     //' compare if accumulated chill is larger than yc
     //' if yes accumulate heat
     if(y[i+1] >= yc){
       //' calculate heat increment and add to z
       z[i+1] = z[i] + P1z(ti, _Tu, _Tb, _Tc);
      
       //' calculate amount of heat needed to achieve flowering (depening on amount of chill)
       z_crit = b1 + (b2 / exp(b3 * x[i+1]));
     }
     
     //' if z_crit calculated (and not zero as default)
     //' check if it is reached, if so return bloom index
     if(z_crit != 0){
       if(z[i+1] >= z_crit){
         // i+2 for Fortran index convention in R
         bloomindex = i+2;
         if(stopatzc) break;
       }
     } 
   }
   if(basic_output) {
     return List::create(Named("bloomindex") = bloomindex);
   }
   return List::create(Named("x") = x, Named("y") = y, Named("z")=z, Named("xs")=xs, Named("bloomindex") = bloomindex);
 }
