% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval_function_endodormancy_meigo.R
\name{eval_function_endodormancy_meigo}
\alias{eval_function_endodormancy_meigo}
\title{Evaluation function for endodormancy}
\usage{
eval_function_endodormancy_meigo(x, SeasonList, dispersion_fun = "calc_cv")
}
\arguments{
\item{x}{model parameters of Dynamic Model, new format: theta_star, theta_c, pie_c, tau, Tf, slope
will get converted to the format: E0, E1, A0, A1, Tf, slope to run the Dynamic Model}

\item{SeasonList}{List of SeasonLists. For each cultivar one entry. Individual seasonlists 
per cultivar contain exactly period of: measured end of endodormancy release to measured
ecodormancy release.}

\item{dispersion_fun}{name of function to measure dispersion of accumulated GDH per
cultivar. Needs to be a function loaded in the environment. Can also be a custom function.
By default "calc_cv", which calculates the coefficient of variation.}
}
\value{
sum of the dispersion function. This will be subject to global optimization.
}
\description{
Takes model parameters for Dynamic Model and list of seasons with daily temperature time series.
}
\details{
The function evaluates a set of parameters and calculates accumulated chill portions (CP) for each of the supplied temperature time series.
Then it returns a dispersion metric measuring how different accumulated CP are among the seasons.
The function is usually used in an optimization problem, where the goal is to estimate
Dynamic Model parameters most suitable to minimize the dispersion of accumulated CP.

The function has very specific requirements on the structure of the Seasonlist.
In fact it is a list of Seasonlists created for several cultivars.
For each cultivar, the seasonlist contains exactly the period from endodormancy 
release to observed ecodormancy release.


Used this function for the Olive study of Efe Deger
}
\examples{
\dontrun{
#          theta_star, theta_c, pie_c, tau, Tf, slope
par <-   c(279,        281,      ?      ?    4,  1.6)

#prepare weather data
weather<-fix_weather(KA_weather[which(KA_weather$Year>2004),])
hourtemps<-stack_hourly_temps(weather, latitude=50.4)
seasonList <- genSeasonList(hourtemps$hourtemps, year = 2006)

#get endodormancy / ecodormancy release date
#--> prepare SeasonList for different cultivars

eval_function_endodormancy_meigo(x = par, par = SeasonList)


}
}
\author{
Lars Caspersen, \email{lars.caspersen@uni-bonn.de}
}
