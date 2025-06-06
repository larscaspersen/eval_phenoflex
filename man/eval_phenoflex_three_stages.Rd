% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eval_phenoflex_three_stages.R
\name{eval_phenoflex_three_stages}
\alias{eval_phenoflex_three_stages}
\title{Evaluation function for three phenological stages at once
 
Takes model parameters for PhenoFlex, however, three individual heat requirements for 
each modelled stage individually.}
\usage{
eval_phenoflex_three_stages(
  x,
  modelfn,
  bloomdates_df,
  SeasonList,
  Tc = 36,
  theta_star = 279,
  na_penalty = 365
)
}
\arguments{
\item{x}{model parameters of PhenoFlex, new format: 
c(yc, zc1, zc2, zc3 s1,Tu, theta_c, pie_c, tau, Tf, Tb, slope)
will get converted to the format: E0, E1, A0, A1, Tf, slope to run the Dynamic Model}

\item{modelfn}{function used within the evaluation function to calculate the actual bloomday, often we use
the 'custom_GDH_wrapper' function for that}

\item{bloomdates_df}{data.frame, with four columns: year, firstbloom, fullbloom and budbreak.
The phenological stages have the julian day of the occurence as values.}

\item{SeasonList}{list of hourly temperatures for the individual phonological seasons. Each element should contain a data.frame
with the columns "Temp" (for the hourly temperature) and "JDay" for the corresponding Julian day. Is usually
generated using \link[chillR]{genSeasonList}}

\item{Tc}{numeric, by default 36. Critical temperature (°C) parameter of the GDH model.}

\item{theta_star}{numeric, by default 279. Optimal temperature (K) for chill accumulation of the
Dynamic Model.}

\item{na_penalty}{numeric, by default 365. Penalty for the phenology
prediction function when it fails return a bloom prediction}
}
\value{
list with two elements. First is called 'f' and contains the residual sum of squares of the model. The 
second is 'g' which is the values of the additional model constraints defined in the function.
If the flag for return_pred set TRUE, then it returns bloom dates
}
\description{
The function returns for one season the timing of the three phenological stages. 
E0 to A1 parameters of the chill submodel are replaced by intermediate parameters. 
Parameters Tc and theta_star are fixed, as they contribute little to
model performance. The function is usually called by an optimization algorithm
during model calibration.
}
\examples{
\dontrun{

ncult <- 3
#          theta_star, theta_c, pie_c, tau, Tf, slope
par <-   c(40, 50, 150, 200, 0.5, 26, 281, 24, 45, 4, 1.6)

#prepare weather data
weather<-fix_weather(KA_weather[which(KA_weather$Year>2004),])
hourtemps<-stack_hourly_temps(weather, latitude=50.4)
seasonList <- genSeasonList(hourtemps$hourtemps, year = 2006)

eval_phenoflex_three_stages(x = par, par = SeasonList)
}
}
\author{
Lars Caspersen, \email{lars.caspersen@uni-bonn.de}
}
