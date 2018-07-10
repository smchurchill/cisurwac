#' Calibrate a loglinear slope
#'
#'@describeIn calibrateSlopeInternal
#'
#'@inheritParams makeTarget
#'@inheritParams makeNormalizedGamma
#'
#'@export
calibrateSlope <- function(count, popn, p_cd, pcc, gamma_relation, lb, ub) {
  target <- makeTarget(count = count, popn = popn)
  mass <- makeNormalizedGamma(p_cd = p_cd, pcc = pcc, gamma_relation = gamma_relation, lb = lb, ub = ub)
  calibrateSlopeInternal(target = target, mass = mass, lb = lb, ub = ub)
}

#' Calibrate a conditional probability density function
#'
#'@inheritParams calibrateSlope
#'
#'@export
calibrateConditionalProbability <- function(count, popn, p_cd, pcc, gamma_relation, lb, ub) {
  target <- makeTarget(count = count, popn = popn)
  mass <- makeNormalizedGamma(p_cd = p_cd, pcc = pcc, gamma_relation = gamma_relation, lb = lb, ub = ub)
  slope <- calibrateSlopeInternal(target = target, mass = mass, lb = lb, ub = ub)
  makeConditionalProbability(slope = slope, lb = lb)
}

#' Calibrate a harm distribution function
#'
#'@inheritParams calibrateSlope
#'
#'@export
calibrateHarmDistribution <- function(count, popn, p_cd, pcc, gamma_relation, lb, ub) {
  target <- makeTarget(count = count, popn = popn)
  rtarget <- 1 / target
  mass <- makeNormalizedGamma(p_cd = p_cd, pcc = pcc, gamma_relation = gamma_relation, lb = lb, ub = ub)
  slope <- calibrateSlopeInternal(target = target, mass = mass, lb = lb, ub = ub)
  risk <- makeConditionalProbability(slope = slope, lb = lb)
  makeHarmDistribution(mass = mass, risk = risk, lb = lb, ub = ub, rtarget = rtarget)
}
