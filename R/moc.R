set.seed(1001)

norm.vec <- function(x) {
  nm.1 <- sum(x)
  return(sapply(x, function (e) { e / nm.1 }))
}


Y <- matrix(c(c(0, 0, 0, 147, 131, 142, 99, 106, 52, 58),
              c(77, 37, 24, 89, 75, 62, 30, 30, 0, 0),
              c(123, 91, 44, 35, 68, 24, 40, 32, 45, 56)),
            nrow=3, byrow=TRUE)
Theta.est <- t(apply(Y, 1, norm.vec))
phi <- c(0.5, 0.5, 0)


## A closure to capture some vars
create.obj.func.f <- function(Theta, comm.idx, phi, y) {
  
  obj.func <- function(x) {
    n <- length(x)
    theta <- x[1:n-1]
    lambda <- x[n]

    Theta.mod <- Theta
    Theta.mod[comm.idx,] <- theta
    mixed.props <- phi %*% Theta.mod
    log.mixed.props <- t(sapply(mixed.props, log))

    penalty <- 3.4^(max(0, sum(theta)))
    result <- (-(log.mixed.props %*% y)) + (penalty * lambda * (sum(theta) - 1))

    return(result)
  }

  return(obj.func)
}


of.f <- create.obj.func.f(Theta.est, 1, phi, Y[1,])

num.dims <- ncol(Y)
result.f <- optim(par=c(Theta.est[1,], 50),
                  fn=of.f,
                  gr=NULL,
                  method="L-BFGS-B",
                  lower=rep(0.0001, num.dims+1),
                  upper=c(rep(1, num.dims), Inf))


create.obj.func.h <- function(Theta, comm.idx, phi, y) {
  
  obj.func <- function(x) {
    n <- length(x)
    theta <- x[1:n-1]
    lambda <- x[n]

    Theta.mod <- Theta
    Theta.mod[comm.idx,] <- theta
    mixed.props <- phi %*% Theta.mod
    
    penalty <- 1000^(max(0, sum(theta)))
    result <- sum((-((phi[comm.idx] / mixed.props) * y) + lambda)^2) + (sum(theta) - 1)^2
    
    return(result)
  }

  return(obj.func)
}

of.h <- create.obj.func.h(Theta.est, 1, phi, Y[1,])
num.dims <- ncol(Y)
result.h <- optim(par=c(Theta.est[1,], 50),
                  fn=of.h,
                  gr=NULL,
                  method="L-BFGS-B",
                  lower=rep(0.0001, num.dims+1),
                  upper=c(rep(1, num.dims), Inf))


create.obj.func.h.2 <- function(Theta, comm.idx, phi, y) {
  
  obj.func <- function(x) {
    n <- length(x)
    theta <- x[1:n]

    Theta.mod <- Theta
    Theta.mod[comm.idx,] <- theta
    mixed.props <- phi %*% Theta.mod

    result <- sqrt(sum((-((phi[comm.idx] / mixed.props) * y))^2))
    return(result)
  }

  return(obj.func)
}

## Like h but with Lagrange term added after finding partials
create.obj.func.h.3 <- function(Theta, comm.idx, phi, y) {
  
  obj.func <- function(x) {
    n <- length(x)
    theta <- x[1:n-1]
    lambda <- x[n]

    Theta.mod <- Theta
    Theta.mod[comm.idx,] <- theta
    mixed.props <- phi %*% Theta.mod
    
    ##penalty <- 1000^(max(0, sum(theta)))
    result <- sum((-((phi[comm.idx] / mixed.props) * y))^2) + lambda * (sum(theta) - 1)
    
    return(result)
  }

  return(obj.func)
}

of.h.3 <- create.obj.func.h.3(Theta.est, 1, phi, Y[1,])
num.dims <- ncol(Y)
result.h.3 <- optim(par=c(Theta.est[1,], 50),
                    fn=of.h.3,
                    gr=NULL,
                    method="L-BFGS-B",
                    lower=c(rep(0.0001, num.dims), 5),
                    upper=c(rep(1, num.dims), Inf))


##Theta.est <- matrix(rep(1/10, 30), nrow=3, byrow=TRUE)

## thetas.community.weights <- matrix(c(norm.vec(c(0, 0, 0, 0.075, 0.05, 0.05, 0.025, 0.025, 0, 0)),
##                                        norm.vec(c(0, 0, 0, 0.025, 0.050, 0.050, 0.075, 0.075, 0.100, 0.100)),
##                                        norm.vec(c(0.3, 0.2, 0.1, 0.05, 0.05, 0, 0, 0, 0, 0))),
##                                        nrow=3, byrow=TRUE)
## Theta.est <- thetas.community.weights
