# An R implementation of Algorithm 1 from BCM17
# Written by Martin S. Copenhaver (www.mit.edu/~mcopen)


FA <- function(S, factors, maxiter.inner = 1000, maxiter.outer = 1000, rho = .01, tol.inner = 1e-9, tol.outer = 1e-5){
  # Input: a covariance matrix S and the desired number of factors to perform factor analysis.
  # Other parameters: optimal algorithmic parameters with defaults
  # Output: decomposition S = T + P + N, where T is positive-semidefinite (PSD) with rank <= factors,
  #         P is diagonal and PSD, and N is PSD (N is the "noise" component)
  
  ### algorithmic parameters:
  # maxiter.* : maximum number of * iterations 
  # rho : scaling parameter in ADMM
  # tol.* : * optimality tolerance
  
  # verify that S is indeed a matrix
  
  if (!(is.matrix(S))) {
    stop(simpleError("Inputted covariance matrix is not in matrix format. Convert using as.matrix( )."));
  }
  
  # verify that S is square
  
  if (dim(S)[1] != dim(S)[2]) {
    stop(simpleError("Inputted covariance matrix is not square."));
  }
  
  
  # verify (approximate symmetry of S)
  
  if (norm(S-t(S), type="F") > 1e-10) {
    stop(simpleError("Inputted covariance matrix is not symmetric."));
  }
  
  # problem parameters
  p = dim(S)[1];
  
  # verify that number of factors is valid (between 0 and p, inclusive)
  
  if (!(as.integer(factors) == factors | factors < 0 | factors > p)) {
    stop(simpleError("Number of factors is not valid."));
  }
  

  # key variables
  phi = matrix(0, nrow = p, ncol = 1);
  W = matrix(1, nrow = p, ncol = p);
  
  # useful constants
  zp = matrix(0, nrow = p, ncol = 1); #matrix of zeros
  dS = diag(S); # diagonal of S
  weig = diag(c(rep(0,factors),rep(1,p-factors))); #weighting used on eigenvalues
  
  # "inner" phi, Lambda, and Nu (for inner iterations)
  lphi = rep(0, p);
  lLambda = matrix(0, nrow = p, ncol = p);
  lNu = matrix(0, nrow = p, ncol = p);
  
  bobj = Inf; # best objective
  
  oerr = Inf; # "outer" error for problem wrt W
  opobj = Inf; # previous (outer) objective
  ocobj = 0; # current (outer) objective
  
  ocount = 0; # outer counter
  while ( (oerr > tol.outer) & (ocount < maxiter.outer) ) {
    ### solve wrt phi
    
    ierr = Inf; # "inner" error for problem wrt Phi (which requires an ADMM)
    pobj = Inf; # previous objective
    cobj = 0; # current objective
    
    dW = diag(W);
    dWs = dW/rho;
    
    count = 0;
    while ( (ierr > tol.inner) & (count < maxiter.inner) ) {
      ## first: update lphi -- closed form
      
      lphi = pmax(zp,dS - diag(lLambda) + dWs - diag(lNu)/rho);
      
      
      ## now: update lLambda -- eigendecomposition (projection on PSD)
      
      e = eigen(S - diag(lphi) - lNu/rho,T); # "T" for symmetric argument
      
      lLambda = e$vectors %*% diag(as.vector(pmax(zp,as.vector(e$values)))) %*% t(e$vectors) # would normally need inverse, but can just do tranpose since unitary.
      
      
      ## then: update lNu -- closed form (gradient update)
      
      lNu = lNu + rho * (lLambda - S + diag(as.vector(lphi)));
      
      
      ## finally: update ierr (error) = norm of (lLambda - (S-lph))
      
      pobj = cobj;
      cobj = sum(dW*lphi);
      relimp = abs(cobj-pobj)/(abs(pobj)+.01)*100;
      ierr = max( norm( lLambda - S + diag(as.vector(lphi)) , type = "F") , relimp );
      # print(c(ierr,count));
      count = count + 1;
      
    }
    
    # update phi to be lphi
    
    phi = lphi;
    
    ### solve wrt W --- requires eigendecomp of S - diag(phi)
    
    e = eigen(S - diag(as.vector(phi)), T );
    W = e$vectors %*% weig %*% t(e$vectors); # would normally need inverse, but can just do tranpose since hermitian (real unitary)
    
    bobj = min(bobj, sum(W*(S-diag(as.vector(phi)))) ); # second term is the inner product of W and S-Phi
    
    opobj = ocobj;
    ocobj = sum(W*(S-diag(as.vector(phi))));
    
    oerr = abs(ocobj-opobj)/(abs(opobj)+.01)*100;
    ocount = ocount + 1;
    
    print(c(ocount,bobj));
  }
  
  e = eigen(S - diag(as.vector(phi)), T );
  Theta = e$vectors %*% diag(c(e$values[1:factors],rep(0,p-factors))) %*% t(e$vectors); # would normally need inverse, but can just do tranpose since hermitian (real unitary)
  
  return(list(Theta=Theta,Phi=as.vector(phi))); # return matrix Theta (T) and Phi (P)
}
