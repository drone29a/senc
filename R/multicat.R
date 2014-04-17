source("moc.R")

set.seed(1001)

Curry <- function(FUN, ...) {
  args <- match.call(expand.dots = FALSE)$...
  args$... <- as.name("...")

  env <- new.env(parent = parent.frame())

  if (is.name(FUN)) {
    fname <- FUN
  } else if (is.character(FUN)) {
    fname <- as.name(FUN)
  } else if (is.function(FUN)){
    fname <- as.name("FUN")
    env$FUN <- FUN
  } else {
    stop("FUN not function or name of function")
  }
  curry_call <- as.call(c(list(fname), args))

  f <- eval(call("function", as.pairlist(alist(... = )), curry_call))
  environment(f) <- env
  f
}


norm.vec <- function(x) {
  nm.1 <- sum(x)
  return(sapply(x, function (e) { e / nm.1 }))
}


## I had started this as 0.7, 0.3 and 0.5, 0.5
thetas <- matrix(c(c(0.7, 0.3),
                   c(0.3, 0.7)),
                 nrow=2, byrow=TRUE)
phi <- c(0.5, 0.5)

## This is the result of using comm. membership to
## calculate the mixed probability that we expect
## is the underlying categorical for this object.
mixed.theta <- phi %*% thetas

## observations, we can think of this as either we use
## one of two community distributions with some probability
## as given by the phi. Or we can think of this is a collapsed
## categorical distribution. Note the latent variable is being
## represented by membership mixture.
## NOTE: Is this intuition correct?
y <- colSums(t(rmultinom(1000, 1, mixed.theta)))

#### Know truth for other communities
## Observe what happens if we estimate when we know the true parameters
## for other communiites we mix with.
## We would usually find a weighted average of observations for all objects
## that are affiliated with the community we are estimating for.
## Since we only have one object, we will use the proportion of the
## observations as our initial estimate.
theta.est <- norm.vec(y)

ll <- function(y, thetas, theta.est.1) {
  theta.est <- c(theta.est.1, 1-theta.est.1)
  thetas.est <- thetas
  thetas.est[1,] <- theta.est
  mixed.probs <- phi %*% thetas.est
  ## NOTE: Is there a problem with using the multinomial likelihood to estimate
  ## for categorical distribution? what happens if we remove multinomial coefficient?
  ## Our current assumption is that the multinomial coefficient is constant and thus can
  ## be dropped for optimization. I.e., the minimum is found at the same point.
  return (-dmultinom(y, prob=mixed.probs, log=TRUE))
}

theta.est.1s <- seq(0.001, 0.999, length.out=100000)
lls <- sapply(theta.est.1s, Curry(ll, y, thetas))

lls.df <- data.frame(theta.est.1=theta.est.1s,
                     log.likelihood=lls)
ggplot(lls.df, aes(x=theta.est.1, y=log.likelihood)) + geom_point()

## The minimum
theta.est.1.min <- theta.est.1s[which.min(lls)]
theta.est.best <- c(theta.est.1.min, 1-theta.est.1.min)

## This min estimate may be a litle off due to randomness, but it comes pretty close.
## E.g., 0.74, 0.26 instead of the actual 0.7, 0.3

#### Both community parameters are estimated.
## Now what happens if we don't know the true parameters for other community
thetas.unk <- matrix(c(theta.est, theta.est), nrow=2, byrow=TRUE)
lls.unk <- sapply(theta.est.1s, Curry(ll, y, thetas.unk))

lls.unk.df <- data.frame(theta.est.1=theta.est.1s,
                         log.likelihood=lls.unk)
ggplot(lls.unk.df, aes(x=theta.est.1, y=log.likelihood)) + geom_point()

## The minimum
theta.est.1.min <- theta.est.1s[which.min(lls.unk)]
theta.est.best <- c(theta.est.1.min, 1-theta.est.1.min)

## NOTE: This is a pathological case! There is only one object shared between two communities.
## In this case, the estimate for each community is the same and since it's based on a single object the likelihood
## is trivially maxmized for that object's features.


#### Instead of using estimate on object, consider using uniform prob on all features
thetas.unk <- matrix(c(c(0.5, 0.5), c(0.6, 0.4)), nrow=2, byrow=TRUE)
lls.unk <- sapply(theta.est.1s, Curry(ll, y, thetas.unk))

lls.unk.df <- data.frame(theta.est.1=theta.est.1s,
                         neg.log.likelihood=lls.unk)

## The minimum
theta.est.1.min <- theta.est.1s[which.min(lls.unk)]
theta.est.best <- c(theta.est.1.min, 1-theta.est.1.min)

ggplot(lls.unk.df, aes(x=theta.est.1, y=neg.log.likelihood)) + geom_point() + geom_vline(xintercept=0.7, linetype="solid") + geom_vline(xintercept=0.486, linetype="dotted") + geom_vline(xintercept=theta.est.best[1], linetype="dashed")



## Now we see clearly that our estimate moves in a direction to account for the difference between the estimated
## parameters for the other communities and their actualy value. This makes it clear that we DO want to estimate
## initial values from the data to avoid turning our community estimation into accounting for the error
## that exists in the mixture of all other communities.

## The conclusion is that we should estimate the community params from the objects and if we only have objects
## that are all members of the same communities, the estimates for all the communities will be the same and we
## will not be able to improve out estimate. That is, our method relies on objects being generated from different
## mixtures of communities.

## Next we should consider multiple objects. One case where the objects are mixtures of the same community but with different mixture weights. A second case where two objects share one community but are each mixed with one other different community. This second case can be broken down further to using both objects together to estimate community params and do a sweep of the ll, as well as looking at each object separately as we have been doing so far. It could be that high dimensions and sparseness of many features will help us when estimating?

#### Two objects with different mixture weights and community clique cores
## Treat each object as core for a community, assume membership is significantly stronger for other
## community
phi.1 <- c(0.75, 0.25)
phi.2 <- c(0.25, 0.75)

y.1 <- colSums(t(rmultinom(1000, 1, phi.1 %*% thetas)))
y.2 <- colSums(t(rmultinom(1000, 1, phi.2 %*% thetas)))

theta.1 <- norm.vec(y.1)
theta.2 <- norm.vec(y.2)
lls.1 <- sapply(theta.est.1s, Curry(ll, y.1, matrix(c(theta.1, theta.2), nrow=2, byrow=TRUE)))

lls.1.df <- data.frame(theta.est.1=theta.est.1s,
                       log.likelihood=lls.1)
ggplot(lls.1.df, aes(x=theta.est.1, y=log.likelihood)) + geom_point()

## The minimum
theta.est.1.min <- theta.est.1s[which.min(lls.1)]
theta.est.best <- c(theta.est.1.min, 1-theta.est.1.min)

lls.2 <- sapply(theta.est.1s, Curry(ll, y.2, matrix(c(theta.2, theta.1), nrow=2, byrow=TRUE)))

lls.2.df <- data.frame(theta.est.1=theta.est.1s,
                       log.likelihood=lls.2)
ggplot(lls.2.df, aes(x=theta.est.1, y=log.likelihood)) + geom_point()

## The minimum
theta.est.1.min <- theta.est.1s[which.min(lls.1)]
theta.est.best <- c(theta.est.1.min, 1-theta.est.1.min)


#### Gradient
thetas <- matrix(c(c(0.7, 0.3),
                   c(0.3, 0.7)),
                 nrow=2, byrow=TRUE)
phi <- c(0.5, 0.5)
y <- colSums(t(rmultinom(1000, 1, phi %*% thetas)))
h <- create.obj.func.h.2(thetas, 1, phi, y)
thetas.sweep <- t(sapply(theta.est.1s, function(x) { c(x, 1-x) }));
mags <- apply(thetas.sweep, 1, h)

mags.df <- data.frame(mag=mags, term.1=thetas.sweep[,1])
ggplot(mags.df, aes(x=term.1, y=mags)) + geom_point()

## Next up, plot the gradient

