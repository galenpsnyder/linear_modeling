mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv") # for logistic

mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/poisson_sim.csv") # for poisson
mydata <- within(mydata, {
  prog <- factor(prog, levels=1:3, labels=c("General", "Academic", 
                                            "Vocational"))
  id <- factor(id)
})

source(file = "/Users/galensnyder/Documents/R scripts/nelder.mead.bounds 101419.R", local = T)

## model parcing module------------------------------------
lm.model.mod <- function(formula, dat, intercept = TRUE){
  ## initializing and formatting function call
  dat <- as.data.frame(dat)
  df <- model.frame(formula, data = dat)
  N <- nrow(df)
  
  namesIV <- names(df)[1]
  namesFE <- names(df)[-1]
  df <- as.matrix(df)
  
  ## formatting response vector
  y <- matrix(df[,1], nrow = N, ncol = 1)
  
  ## formatting predictor matrix
  x <- matrix(df[, -1], nrow = N)
  if(ncol(x) == 0){
    hessian = FALSE; cat("warning: hessian unreliable for 1 parameter model\n")
    empty <- TRUE
    p <- 1
    x <- matrix(1, nrow = N)
    namesFE <- "(intercept)"
  } else {
    hessian = TRUE
    empty <- FALSE
    if(intercept == TRUE){
      x <- cbind(1, x)
      p <- ncol(x)
      namesFE <- c("(intercept)", namesFE)
    }
  }
  
  mod <- list(
    x = x, y = y, N = N, p = p, intercept = intercept, empty = empty, namesFE = namesFE, namesIV = namesIV, hessian = hessian
  )
  mod
}


## optimization module-----------------------------
flexlm.fit <- function(theta, mod = NULL, x = NULL, y = NULL, intercept = TRUE, family = "gaussian", tol = .Machine$double.eps){
  if(!is.null(mod)){
    N <- mod$N
    y <- mod$y
    x <- mod$x
    p <- mod$p
  } else {
    y <- y
    x <- x
    N <- length(y)
    p <- ncol(x)
  }
  
  
  if(family == "gaussian"){
    beta  <- theta
    y_hat <- x%*%beta
    rss <- sum((y - y_hat)^2)
    fit <- -2*(0.5*(-N*(log(2*pi) + 1 - log(N) + log(rss)))) 
  }
  
  if(family == "binomial"){
    beta  <- theta
    y_hat <- x%*%beta
    fit <- -2 * sum(y*y_hat - log(1 + exp(y_hat)))
  }
  
  if(family == "poisson"){
    beta  <- theta
    y_hat <- x%*%beta
    fit <- -2 * sum(y*y_hat - exp(y_hat))
  }
  
  if(family == "multinom"){
    if(!is.factor(y)) y <- as.factor(y)
    k     <- length(levels(y))
    beta  <- matrix(theta, p, k-1)
    new_y <- matrix(0, nrow = length(y), ncol = k)
    new_y <- cbind(y, new_y)
    for(i in 1:k){
      for(j in 1:N){
        ifelse(new_y[j, 1] == i, new_y[j, i+1] <- 1, new_y[j, i+1] <- 0)
      }
    }
    new_y <- new_y[, -c(1:2)]
    y_hat <- x%*%beta
    fit <- -2*sum(rowSums(new_y*y_hat)-log(1 + rowSums(exp(y_hat))))
  }
  fit
}


## output module-------------------------------
lm.output.mod <- function(fit, mod, x = NULL, y = NULL, family = NULL, intercept = TRUE, hessian = TRUE, return.hessian = FALSE){
  beta    <- b <- fit$parm
  dev     <- fit$value
  LL      <- -0.5*dev
  iter    <- fit$iterations
  x       <- mod$x
  y       <- mod$y
  N       <- mod$N
  p       <- mod$p
  namesFE <- mod$namesFE
  empty   <- mod$empty
  pred    <- x %*% beta
  aic     <- dev + 2*p
  bic     <- dev + p*log(N)
  
  if(hessian == TRUE){
    hess  <- 0.5*fit$hessian
    var.b <- sqrt(diag(solve(hess)))
  } else {
    hess  <- NULL
    var.b <- rep(NA, p)
  }
  
  if(family == "gaussian"){
    ifelse(!is.null(hess), t.val <- beta / var.b, t.val <- rep(NA, p))
    rss   <- sum((y - pred)^2)
    ifelse(intercept == T, pf <- p-1, pf <- p)
    ifelse(intercept == T, diff  <- pred - mean(y), diff <- pred)
    ssr   <- sum(diff^2)
    sigma <- rss / (N-p)
    f.val <- (ssr / (pf)) / sigma
    r2    <- 1 - (rss / (rss + ssr))
    adj.r <- 1 - (1-r2)*((N-1)/(N-p))
    
    beta <- cbind(beta, var.b, t.val)
    beta <- data.frame(beta, row.names = namesFE)
    names(beta) <- c("Estimate", "Std. error", "t value")
  }
  
  if(family == "binomial"){
    ifelse(!is.null(hess), z.val <- beta / var.b, z.val <- rep(NA, p))
    beta <- cbind(beta, var.b, z.val)
    beta <- data.frame(beta, row.names = namesFE)
    names(beta) <- c("Estimate", "Std. error", "z value")
  }
  
  if(family == "poisson"){
    if(!is.null(hess)){
      hess <- hess/2
      var.b   <- sqrt(diag(solve(hess)))
      z.val <- beta / var.b
    } else {
      z.val <- rep(NA, p)
    }
    
    LL   <- LL - sum(log(factorial(y)))
    dev  <- -2*LL
    aic  <- dev + 2*p
    bic  <- dev + p*log(N)
    beta <- cbind(beta, var.b, z.val)
    beta <- data.frame(beta, row.names = namesFE)
    names(beta) <- c("Estimate", "Std. error", "z value")
  }
  
  
  cat("General linear model fit by ML estimation\n")
  cat("Statistical model\n")
  cat("  Family:", family, "\n") 
  cat("  Number of observations:", N, "\n")
  cat("  Number of parameters:", p, "\n")
  cat("\n")
  cat("Iteration process:\n")
  cat("  Model converged normally after", iter, "iterations\n")
  # cat("  Time elapsed:", unname(etm[3]), "seconds\n")
  cat("\n")
  cat("Goodness-of-fit statistics\n")
  cat("  Log-likelihood:", LL, "\n")
  cat("  AIC:", aic, "\n")
  cat("  BIC:", bic, "\n")
  cat("\n")
  cat("Coefficients\n")
  print(beta)
  cat("\n")
  if(family == "gaussian"){
    cat("Residual variance:", sigma, "on", N-p, "degrees of freedom\n")
    if(empty == FALSE) cat("R-squared:", r2, "\n")
    if(empty == FALSE) cat("Adjusted R-squared:", adj.r, "\n")
    if(empty == FALSE) cat("F-statistic:", f.val, "on", pf, "and", N-p, "degrees of freedom\n") 
  }

  if(family == "binomial" || family == "poisson"){
    cat("Null deviance:", mod$null - mod$sat, "on", N-1,  "degrees of freedom\n")
    cat("Residual deviance:", dev - mod$sat, "on", N-p, "degrees of freedom\n")
    cat("-----\n")
  }
  
  out <- list(
    fitted     = pred,
    parameters = b,
    loglik     = LL
  )
  invisible(out)
}



flexlm <- function(formula = NULL, dat = NULL, intercept = TRUE, family = "gaussian", return.hessian = FALSE){
  if(is.null(formula)) stop("No formula specified!")
  
  mod <- lm.model.mod(formula = formula, dat = dat, intercept = intercept)
  
  if(family == "binomial" && any(mod$y %in% c(0, 1) == FALSE)){
    stop(gettext("y values must be coded 0, 1"), call. = F)
  }
  
  if(min(mod$y) > 0){
    check.dist<- TRUE
    fam.len   <- max(length(mod$y), 100)
    fam.check <- quantile(mod$y, seq(0.01, .99, length.out = fam.len))
    fam.refn  <- qnorm(p = seq(0.01, .99, length.out = fam.len), mean = mean(mod$y), sqrt(var(mod$y)))
    fam.refp  <- qpois(p = seq(0.01, .99, length.out = fam.len), mean(mod$y))
    
    fam.refn  <- sum((fam.check - fam.refn)^2)
    fam.refp  <- sum((fam.check - fam.refp)^2)
  } else {
    check.dist <- FALSE
  }
  
  
  
  if(family == "binomial" || family == "poisson" || family == "multinom"){
    null.fit <- nelder.mead(1, flexlm.fit, tol = .Machine$double.eps, 
                            y = mod$y, x = matrix(1, nrow = mod$N), family = family)
    mod$null <- null.fit$value
    mod$sat <- 0
    
    if(family == "poisson"){
      mod$null <- mod$null + 2*sum(log(factorial(mod$y)))
      mod$sat <- -2*(sum(mod$y[mod$y != 0] * log(mod$y[mod$y != 0])) - sum(mod$y[mod$y != 0]) -
        sum(log(factorial(mod$y[mod$y != 0]))))
    }
  }
  theta.start <- runif(mod$p, 0.1, .9)
  
  if(check.dist && family == "gaussian" && fam.refp < fam.refn){
    hist(mod$y, main = "Histogram of Outcome Variable", xlab = mod$namesIV)
    warning(gettext("Family == 'poisson' may be more appropriate for outcome variable..."), call. = F, immediate. = T)
    user.response <- readline("Would you like to estimate alternative model for comparison? (yes/no): ")
    while(user.response != "yes" && user.response != "no"){
      user.response <- readline("Please type yes or no to proceed: ")
    }
    if(user.response == "yes"){
      cat("-----ALTERNATIVE MODEL-----\n")
      robust.check.fit <- nelder.mead(theta.start, flexlm.fit, hessian = mod$hessian, tol = .Machine$double.eps, mod = mod, intercept = intercept, family = "poisson")
      lm.output.mod(fit = robust.check.fit, mod = mod, intercept = intercept, family = "poisson", hessian = mod$hessian)
      cat("\n")
      cat("-----USER SPECIFIED MODEL-----\n")
    }
  }
  
  if(family == "poisson" && fam.refn < fam.refp){
    hist(mod$y, main = "Histogram of Outcome Variable", xlab = mod$namesIV)
    warning(gettext("Family == 'gaussian' may be more appropriate for outcome variable..."), call. = F, immediate. = T)
    user.response <- readline("Would you like to estimate alternative model for comparison? (yes/no): ")
    while(user.response != "yes" && user.response != "no"){
      user.response <- readline("Please type yes or no to proceed: ")
    }
    if(user.response == "yes"){
      cat("-----ALTERNATIVE MODEL-----\n")
      robust.check.fit <- nelder.mead(theta.start, flexlm.fit, hessian = mod$hessian, tol = .Machine$double.eps, mod = mod, intercept = intercept, family = "gaussian")
      lm.output.mod(fit = robust.check.fit, mod = mod, intercept = intercept, family = "gaussian", hessian = mod$hessian) 
      cat("\n")
      cat("-----USER SPECIFIED MODEL-----\n")
    }
  }
  
  fit <- nelder.mead(theta.start, flexlm.fit, hessian = mod$hessian, tol = .Machine$double.eps, mod = mod, intercept = intercept, family = family)
  
  lm.output.mod(fit = fit, mod = mod, intercept = intercept, family = family, hessian = mod$hessian)
}

flexlm(formula = admit~gpa, dat = mydata, family = "binomial")
flexlm(formula = num_awards ~ math, dat = mydata, family = "poisson")
flexlm(formula = gre~gpa, dat = mydata)

