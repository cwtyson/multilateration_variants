
invest(object = rssi_dist_mdl, 
       y0 = -10,
       mean.response = FALSE,
       level=0.95,
       data = mdl_data)

.data  <- if (!missing(data)) data else eval(object$call$data, 
                                             envir = parent.frame())
x0.name <- intersect(all.vars(stats::formula(object)[[3]]), colnames(.data)) 
yo = -73
m <- length(y0)  # number of unknowns 
eta <- mean(y0)  # mean response
n <- length(stats::resid(object))  # sample size
p <- length(stats::coef(object))  # number of parameters
var.pooled <- stats::sigma(object)^2  # residual variance
lower <- min(.data[, x0.name])  # lower limit default
upper <- max(.data[, x0.name])  # upper limit default
x = .data
tol = .Machine$double.eps^0.25 
maxiter = 1000

  
computeInverseEstimate <- function(object, x0.name, eta, lower, upper, 
                                       extendInt, tol, maxiter) {

  # Calculate point estimate by inverting fitted model
  res <- try(stats::uniroot(function(x) {
    stats::predict(object, newdata = makeData(x, x0.name)) - eta
  }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
  silent = TRUE)
  
  # Provide (informative) error message if point estimate is not found
  if (inherits(res, "try-error")) {
    stop(paste("Point estimate not found in the search interval (", lower, 
               ", ", upper, "). ", 
               "Try tweaking the values of lower and upper. ",
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  } else {
    res
  }
  
}

makeData <- function(x, label) {
  stats::setNames(data.frame(x), label)
}  

# Calculate point estimate by inverting fitted model
x0.est <- computeInverseEstimate(object, x0.name = x0.name, eta = eta, 
                                 lower = lower, upper = upper, 
                                 extendInt = extendInt, tol = tol,
                                 maxiter = maxiter)
