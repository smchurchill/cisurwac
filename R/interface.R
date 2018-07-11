#' Calibrate a loglinear slope
#'
#'@describeIn calibrateSlopeInternal
#'
#'@inheritParams makeTarget
#'@inheritParams makeNormalizedGamma
#'
#'@export
calibrateSlope <- function(count, popn, p_cd, pcc, gamma_relation, lb, ub, threshold) {
  target <- makeTarget(count = count, popn = popn, p_cd = p_cd)
  mass <- makeNormalizedGamma(p_cd = p_cd, pcc = pcc, gamma_relation = gamma_relation, lb = lb, ub = ub)
  calibrateSlopeInternal(target = target, mass = mass, threshold = threshold, ub = ub)
}

#' Calibrate a cumulative conditional probability function
#'
#'@inheritParams calibrateSlope
#'
#'@export
calibrateConditionalProbability <- function(count, popn, p_cd, pcc, gamma_relation, lb, ub, threshold) {
  target <- makeTarget(count = count, popn = popn, p_cd = p_cd)
  mass <- makeNormalizedGamma(p_cd = p_cd, pcc = pcc, gamma_relation = gamma_relation, lb = lb, ub = ub)
  slope <- calibrateSlopeInternal(target = target, mass = mass, threshold = threshold, ub = ub)
  makeConditionalProbability(slope = slope, threshold = threshold)
}

#' Calibrate a harm distribution function
#'
#'@inheritParams calibrateSlope
#'
#'@export
calibrateHarmDistribution <- function(count, popn, p_cd, pcc, gamma_relation, lb, ub, threshold) {
  target <- makeTarget(count = count, popn = popn, p_cd = p_cd)
  mass <- makeNormalizedGamma(p_cd = p_cd, pcc = pcc, gamma_relation = gamma_relation, lb = lb, ub = ub)
  slope <- calibrateSlopeInternal(target = target, mass = mass, threshold = threshold, ub = ub)
  risk <- makeConditionalProbability(slope = slope, threshold = threshold)
  makeHarmDistribution(mass = mass, risk = risk, lb = lb, ub = ub, target = target)
}
