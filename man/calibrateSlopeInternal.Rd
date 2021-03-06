% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interface.R, R/internal.R
\name{calibrateSlope}
\alias{calibrateSlope}
\alias{calibrateSlopeInternal}
\title{Calibrate a loglinear slope}
\usage{
calibrateSlope(count, popn, p_cd, pcc, gamma_relation, lb, ub, threshold)

calibrateSlopeInternal(target, mass, threshold, ub)
}
\arguments{
\item{count}{Total number of events among the given population}

\item{popn}{Total number of persons in the given population}

\item{p_cd}{Proportion of given population that are current drinkers}

\item{pcc}{daily per capita consumption in grams}

\item{gamma_relation}{relationship between mean and standard deviation of the
gamma distribution that governs alcohol consumption.  By Kehoe, use 1.171 for
male populations and 1.258 for female populations.}

\item{lb}{lower bound of consumption}

\item{ub}{upper bound of consumption}

\item{threshold}{lower bound of consumption at which condition occurs}

\item{target}{Observed incidence to calibrate against}

\item{mass}{population exposure mass function to calibrate against}
}
\value{
slope of loglinear conditional probability mass function for risk as a
result of exposure
}
\description{
Calibrate a loglinear slope

Given a list that contains the true count among a given cohort of
a wholly alocohol-attributable condition, the number of drinkers in that
cohort, the conditions's IM code, and the gamma function detailing consumption
among the cohort, this function calibrates a slope for a loglinear function
that estimates the conditional probablity of developing the given condition.
The function should be interpreted as conditional probability mass.  I.e., the
integral of the function over the relevant range is the probability that a
drinker will be afflicted by an event of the given condition over the period
of time that the given Count variable was collected.

Uses a nonlinear optimizer (COBYLA) to find a loglinear slope for the
function f(x) = max(1, exp(k(x-t))) that mimizes the difference between
integral(N_GAMMA * (f-1), 0.03, UB) and yearly prevalence (count/drinkers).

The goal is to produce a continuous analogue to the relative risk curve for
conditions that are wholly attributable to alcohol. The assumption is made
that such a condition has a loglinear thresholded (i.e. f(x)=1 for x<t)
conditional probablity function on the interval of concern (0.03 to UB grams
of ethanol/day, averaged over 1yr).

This conditional probability is used to portion a 1.00 AAF_TOTAL among the
drinking population.
}
\section{Functions}{
\itemize{
\item \code{calibrateSlope}: 
}}

