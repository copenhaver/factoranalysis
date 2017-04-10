# A demo of R implementation of Algorithm 1 from BCM17
# (as contained in code.R)
# Written by Martin S. Copenhaver (www.mit.edu/~mcopen)


######################

# generate random Theta and Phi (class A_1 in BCM17)

p = 20
r = 2

# create matrices THETA and PHI

L = matrix(rnorm(p*r), ncol = r)
THETA = L %*% t(L)

PHI = runif(p)
PHI = (sum(diag(THETA))/sum(PHI))*PHI # normalized so that equal proportion of common and individual variances

# covariance matrix S is sum of THETA and PHI

S = THETA + diag(PHI)


####  Now perform factor analysis

res = FA(S,2) # FA function from code.R

# See if you recover true Theta and true Phi

norm(THETA - res$Theta, type="F") # in Frobenius norm

norm(PHI - as.matrix(res$Phi), type="F")
