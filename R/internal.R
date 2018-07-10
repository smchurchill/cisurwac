#' Make a target incidence for loglinear slope calibration
#'
#'@param count Total number of events among the given population
#'@param popn Total number of persons in the given population

makeTarget <- function(count, popn) {
  if(count <= 0) stop("Count must be strictly positive.")
  if(popn <= 0) stop("Population must be strictly positive.")
  count / (p_cd * popn)
}


#' Factory for a gamma density normalized over the given interval to the given
#' drinking population
#'
#'@param p_cd Proportion of given population that are current drinkers
#'@param pcc daily per capita consumption in grams
#'@param gamma_relation relationship between mean and standard deviation of the
#'gamma distribution that governs alcohol consumption.  By Kehoe, use 1.171 for
#'male populations and 1.258 for female populations.
#'@param lb lower bound of consumption
#'@param lb upper bound of consumption

makeNormalizedGamma <- function(p_cd, pcc, gamma_relation, lb, ub) {
  if(p_cd <= 0 || p_cd >= 1) stop("Drinking proportion must be between 0 and 1.")
  if(pcc <= 0) stop("Daily consumption must be strictly positive.")
  if(lb < 0 || ub <= lb) stop("Lower bound of consumption must be nonnegative and less than the upper bound of consumption.")
  shape <- 1 / gamma_relation / gamma_relation
  scale <- gamma_relation * gamma_relation * pcc
  nc <- pgamma(q = ub, shape = shape, scale = scale) -
    pgamma(q = lb, shape = shape, scale = scale)
  df = p_cd / nc

  function(x) df * dgamma(x, shape = shape, scale = scale)
}

#' Factory for loglinear conditional probability functions
#'
#'

makeConditionalProbability <- function(slope, lb) {
  function(x) exp(pmax(0, slope*(x-lb)))-1
}

#' Factory for a harm distribution function
#'
#'

makeHarmDistribution <- function(mass, risk, lb, ub, rtarget) {
  function(x) rtarget * makeIntegrator(mass %prod% risk, lb, ub)(x)
}

#' Factory for integrators
#' @export
makeIntegrator <- function(f, lb, ub) {
  integrate_up_to <- function(to) {
    if(to <= lb) to = lb
    if(to >= ub) to = ub
    integrate(f = f, lower = lb, upper = to)$value
  }

  function(x) {
    vapply(x, integrate_up_to, 0)
  }
}

#' Factory for Pointwise Function Products
#'@description factory that produces the product of a pair of functions, where
#'the product used is pointwise multiplication
#'
#'@param f,g function that takes a single argument and produces a value that is
#'a valid argument for the `*` function
#'
#' @export
makeProduct <- function(f, g) {
  function(x) f(x) * g(x)
}

#' Binary operator for Pointwise Function Products
#'
#'@description binary operator for product_factory
#'
#'@describeIn makeProduct
#'
#'@inheritParams makeProduct
#' @export
`%prod%` <- function(f,g) makeProduct(f,g)



#' Calibrate a loglinear Slope
#'
#'@description
#'Given a list that contains the true count among a given cohort of
#'a wholly alocohol-attributable condition, the number of drinkers in that
#'cohort, the conditions's IM code, and the gamma function detailing consumption
#'among the cohort, this function calibrates a slope for a loglinear function
#'that estimates the conditional probablity of developing the given condition.
#'The function should be interpreted as conditional probability mass.  I.e., the
#'integral of the function over the relevant range is the probability that a
#'drinker will be afflicted by an event of the given condition over the period
#'of time that the given Count variable was collected.
#'
#'Uses a nonlinear optimizer (COBYLA) to find a loglinear slope for the
#'function f(x) = max(1, exp(k(x-t))) that mimizes the difference between
#'integral(N_GAMMA * (f-1), 0.03, UB) and yearly prevalence (count/drinkers).
#'
#'The goal is to produce a continuous analogue to the relative risk curve for
#'conditions that are wholly attributable to alcohol. The assumption is made
#'that such a condition has a loglinear thresholded (i.e. f(x)=1 for x<t)
#'conditional probablity function on the interval of concern (0.03 to UB grams
#'of ethanol/day, averaged over 1yr).
#'
#'This conditional probability is used to portion a 1.00 AAF_TOTAL among the
#'drinking population.
#'
#'
#'@param target Observed incidence to calibrate against
#'@param mass population exposure mass function to calibrate against
#'@param lb lower bound of consumption at which condition occurs
#'@param ub upper bound of consumption
#'
#'@return slope of loglinear conditional probability mass function for risk as a
#'result of exposure
#'
#' @export
calibrateSlopeInternal <- function(target, mass, lb, ub) {
  if(is.na(target) | target <= 0) return(0)

  integrand <- function(k) function(x) mass(x) * (exp(pmax(0, k*(x-lb)))-1)
  estimate <- function(k) integrate(integrand(k), lb, ub)$value
  error <- function(k) abs(estimate(k) - target)

  nloptr::nloptr(
    x0 = 0.01,
    eval_f = error,
    lb = 0,
    ub = 1,
    opts = list(
      "algorithm" = "NLOPT_LN_COBYLA",
      "xtol_rel" = 1.0e-20
    )
  )$solution
}
