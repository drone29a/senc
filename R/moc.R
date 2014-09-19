set.seed(1001)

require(MCMCpack)

ddirmult <- function(x, alpha) {
    (gamma(A) / gamma(N + A)) 
}

Y <- matrix(c(c(0, 0, 0, 147, 131, 142, 99, 106, 52, 58),
              c(77, 37, 24, 89, 75, 62, 30, 30, 0, 0),
              c(123, 91, 44, 35, 68, 24, 40, 32, 45, 56)),
            nrow=3, byrow=TRUE)
Theta.est <- t(apply(Y, 1, NormVec))
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


## Obj fun for finding mixture of commuities for observations
create.obj.func.mix <- function(Theta, obs) {
  
  obj.func <- function(x) {
    n <- length(x)
    phi <- x[1:n-1]
    lambda <- x[n]

    
    ##penalty <- 1000^(max(0, sum(theta)))
    result <- sum(((phi %*% Theta) - obs)^2) + (lambda * (sum(Theta) - 1))
    
    return(result)
  }

  return(obj.func)
}

of.mix <- create.obj.func.mix(Theta.est, NormVec(Y[2,]))
num.dims <- ncol(Y)
num.comms <- nrow(Theta.est)
result.mix <- optim(par=c(0.5, 0.499, 0.001, 0.01),
                    fn=of.mix,
                    gr=NULL,
                    method="L-BFGS-B",
                    lower=c(rep(0.0001, num.comms), 0.001),
                    upper=c(rep(1, num.comms), Inf))



##Theta.est <- matrix(rep(1/10, 30), nrow=3, byrow=TRUE)

## thetas.community.weights <- matrix(c(norm.vec(c(0, 0, 0, 0.075, 0.05, 0.05, 0.025, 0.025, 0, 0)),
##                                        norm.vec(c(0, 0, 0, 0.025, 0.050, 0.050, 0.075, 0.075, 0.100, 0.100)),
##                                        norm.vec(c(0.3, 0.2, 0.1, 0.05, 0.05, 0, 0, 0, 0, 0))),
##                                        nrow=3, byrow=TRUE)
## Theta.est <- thetas.community.weights

## Generate some documents and see if we can estimate mixtures for docs.

Y <- matrix(c(c(1, 1, 1, 147, 131, 142, 99, 106, 52, 58),
              c(77, 37, 24, 89, 75, 62, 30, 30, 1, 1),
              c(23, 91, 44, 35, 68, 24, 40, 232, 145, 356)),
            nrow=3, byrow=TRUE)
Theta.est <- t(apply(Y, 1, NormVec))

mixes <- c(rep(c(0.7, 0.2999, 0.0001), 5),
           rep(c(0.2, 0.7999, 0.0001), 5),
           rep(c(0.0001, 0.3499, 0.65), 5))
mixes <- matrix(mixes, nrow=length(mixes)/nrow(Y), byrow=TRUE)

docs.params <- mixes %*% Theta.est
docs <- matrix(apply(docs.params, 1, Curry(rmultinom, 1, 100)),
               nrow=nrow(mixes), byrow=TRUE)

p.d.z <- docs %*% t(log(Theta.est))
p.d.z <- exp(p.d.z)

p.d.z <- apply(p.d.z, 2, NormVec)
p.d.z <- t(apply(p.d.z, 1, NormVec))


rows <- lapply(seq_len(nrow(Theta.est)), function(i) Theta.est[i,])

## P(d | Z)
OldDocProbs <- function(d, Theta) {
    num.topic <- nrow(Theta)
    NormVec(apply(Theta, 1, Curry(dmultinom, d, NULL)))
}

## Weight method
TermWeight <- function(term.probs, term.count) {
    if (term.count == 0) {
        return(rep(0, length(term.probs)))
    }

    NormVec(term.probs) * term.count
}

## P(Z | d)
CommWeight <- function(Theta, d) {
    Theta.cols <- split(Theta, col(Theta))
    rowSums(mapply(TermWeight, Theta.cols, d))
}

## P(Z | D) based on prob-weighted term counts
DocWeight <- function(Theta, docs) {
    num.docs <- length(docs)
    counts <- t(apply(docs, 1, Curry(CommWeight, Theta)))
    t(apply(counts, 1, NormVec))
}

## Best method
TermBest <- function(term.probs, term.count) {
    result <- rep(0, length(term.probs))
    result[which.max(term.probs)] <-  term.count
    result
}

CommBest <- function(Theta, d) {
    Theta.cols <- split(Theta, col(Theta))
    rowSums(mapply(TermBest, Theta.cols, d))
}

DocBest <- function(Theta, docs) {
    num.docs <- length(docs)
    counts <- t(apply(docs, 1, Curry(CommBest, Theta)))
    t(apply(counts, 1, NormVec))    
}

## Data as prior, empirical Bayes
CommEmp <- function(Theta.i, d) {
    dmultinom(d, prob=Theta.i) * ddirichlet(Theta.i, d)
}

MixEmp <- function(Theta, d) {
    NormVec(apply(Theta, 1, function (x) { CommEmp(x, d) }))
}

DocEmp <- function(Theta, docs) {
    t(apply(docs, 1, Curry(MixEmp, Theta)))
}

## Weight method with P(Z|\alpha)
TermWt <- function(term.probs, term.count) {
    if (term.count == 0) {
        return(rep(0, length(term.probs)))
    }

    NormVec(term.probs) * term.count
}

## P(Z | d)
CommWt <- function(Theta, d) {
    Theta.cols <- split(Theta, col(Theta))
    rowSums(mapply(TermWeight, Theta.cols, d))
}

## P(Z | D) based on prob-weighted term counts
DocWt <- function(Theta, docs) {
    num.docs <- length(docs)
    counts <- t(apply(docs, 1, Curry(CommWeight, Theta)))
    t(apply(counts, 1, NormVec))
}


## New idea!
## What if we "count" the number of times a topic appears in a document by
## summing the probabilities of each word given each document?
## Then this "count" is used to update the Dirichlet prior and the mixture
## used for the next iteration is the expected value of the Dirichlet distribution.
## Note that this is different than Bayesian inference which relies on credible
## intervals. We are finding a point estimate. We don't have to do this, btw.
## We could use a Dirichlet-multinomial in the other step, but need to determine
## how this changes the equation for analytically solving.

## The above seems to always lead to sampling. E.g., we have a Dirichlet distribution
## and to find the spread we'd sample from it. And we use those samples in our inference?
## What if instead we update this Dirichlet distribution like before, estimate the topic
## mixtures as we are doing right now, but then also include the probability of the
## direct estimate mixture according to the Dirichlet distribution?

## Another interesting point. That we are not taking advantage of the dependencies
## across dimensions when calculating the probability of observed individual counts.
## We might instead want to use a multinomial distribution instead of P(w|z) since
## it will account for how well a community explains then entire count. But then
## how can we normalize this? Or rather, how do we compare the probabilities of the different
## distributions so that when we normalize them one of them doesn't completely dominate
## because they're each measuring different things?

ggplot(data.frame(dim.id=factor(c(sapply(1:10, Curry(rep, times=3)))),
                  comm.id=factor(rep(1:3, 10)),
                  prob=c(Theta.est)),
       aes(x=dim.id, y=prob, group=comm.id, fill=comm.id)) + geom_bar(stat="identity") + facet_grid(~ comm.id)


#### Quick and dirty data generation for MoC tests

comm.a.params <- c(0.7, 0.2, 0.1, 0.0, 0.0)
comm.b.params <- c(0.7, 0.0, 0.3, 0.0, 0.0)
comm.c.params <- c(0.3, 0.0, 0.2, 0.2, 0.3)
comm.d.params <- c(0.4, 0.5, 0.1, 0.0, 0.0)
comm.e.params <- c(0.1, 0.1, 0.2, 0.4, 0.2)

comm.params <- matrix(c(comm.a.params,
                        comm.b.params,
                        comm.c.params,
                        comm.d.params,
                        comm.e.params),
                      nrow=5, byrow=TRUE)

num.objs.per.comm <- c(60, 30, 10, 5, 3)

num.obs.per.obj <- 100

## 1s and 0s, since we aren't explicitly modeling individual membership
## We could also create random membership vectors and then randomly select some which belong
## to each community.

## Various mixing but all with equal concentration
obj.memb.priors.1 <- matrix(c(1.0, 1.0, 0.0, 0.0, 0.0,
                              1.0, 1.0, 1.0, 0.0, 0.0,
                              0.0, 1.0, 1.0, 0.0, 0.0,
                              1.0, 0.0, 0.0, 1.0, 1.0,
                              0.0, 1.0, 0.0, 1.0, 1.0),
                            nrow=5, byrow=TRUE)

## Various mixing pattern with higher concentration of selected community
obj.memb.priors.2 <- matrix(c(10.0, 1.0, 0.0, 0.0, 0.0,
                              1.0, 10.0, 1.0, 0.0, 0.0,
                              0.0, 1.0, 10.0, 0.0, 0.0,
                              1.0, 0.0, 0.0, 10.0, 1.0,
                              0.0, 1.0, 0.0, 1.0, 10.0),
                            nrow=5, byrow=TRUE)

DrawMemb <- function(prior, n) {
    return(rdirichlet(n=n, alpha=prior))
}

DrawFeatVals <- function(comm.params, obj.memb, num.obs) {
    mix.param <- obj.memb %*% comm.params
    return(rowSums(rmultinom(n=num.obs, size=1, prob=mix.param)))
}

GenMembs <- function(priors, num.objs) {
    membs <- list()
    for (i in 1:length(num.objs)) {
        membs[[i]] <- DrawMemb(priors[i,], num.objs[i])
    }
    return(membs)
}

GenMembs <- function(priors, num.objs) {
    membs <- list()
    for (i in 1:length(num.objs)) {
        membs[[i]] <- DrawMemb(priors[i,], num.objs[i])
    }
    return(membs)
}

GenMembsFlat <- function(priors, num.objs) {
    membs <- matrix(nrow=
    for (i in 1:length(num.objs)) {
        membs <- append(membs, unlist(DrawMemb(priors[i,], num.objs[i])))
    }
    return(unlist(membs))
}

GenFeatVals <- function(comm.params, obj.membs, num.obs) {
    num.comms <- length(obj.membs)
    feat.vals <- list() 
    for (i in 1:num.comms) {
        feat.vals[[i]] <- list()
        for (j in 1:(length(obj.membs[[i]]) / num.comms)) {
            feat.vals[[i]][[j]] <- DrawFeatVals(comm.params, obj.membs[[i]][j,], num.obs)
        }
    }
    return(feat.vals)
}

GenFeatValsFlat <- function(comm.params, obj.membs, num.obs) {
    num.comms <- length(obj.membs)
    feat.vals <- list() 
    for (i in 1:num.comms) {
        for (j in 1:(length(obj.membs[[i]]) / num.comms)) {
            feat.vals <- append(feat.vals, DrawFeatVals(comm.params, obj.membs[[i]][j,], num.obs))
        }
    }
    return(matrix(feat.vals, ncol=5, byrow=TRUE))
}


GenData <- function(comm.params, obj.memb.priors, num.objs.per.comm) {
    membs <- GenMembs(obj.memb.priors, num.objs.per.comm)
    ##feat.vals <- GenFeatVals(comm.params, membs, num.obs.per.obj)
    feat.vals <- GenFeatValsFlat(comm.params, membs, num.obs.per.obj)
    return(list(membs=membs,
                feat.vals=feat.vals))
}

FormatMembs <- function(obj.membs) {
    data <- apply(obj.membs,
                  1,
                  function(nz.ind) {
                      sprintf("%s %s %s", nz.ind[1], nz.ind[2], m[nz.ind[1], nz.ind[2]])
                  })
}

FormatFeatVals <- function(feat.vals) {
    data <- apply(feat.vals,
                  1,
                  function(row) {
                      print(sprintf("%s", row))
                  })
    return(data)
}

### Data was generated with:
    ## set.seed(1001)
    ## > GenData <- function(comm.params, obj.memb.priors, num.objs.per.comm) {
    ## +     membs <- GenMembs(obj.memb.priors, num.objs.per.comm)
    ## +     ##feat.vals <- GenFeatVals(comm.params, membs, num.obs.per.obj)
    ## +     feat.vals <- GenFeatValsFlat(comm.params, membs, num.obs.per.obj)
    ## +     return(list(membs=membs,
    ## +                 feat.vals=feat.vals))
    ## + }
