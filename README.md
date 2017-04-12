### Background

Factor analysis ("FA") is a classical technique in multivariate statistics enjoying a rich history dating back over 100 years. In its simplest form, FA focuses on decomposing a given covariance matrix `S` into two respective components: `S=T+P`, where `T` is a (low-rank) positive semidefinite ("PSD") matrix (representing the common variances) and `P` is a diagonal matrix with non-negative diagonal entries (representing the individual variances). In reality, this *noiseless* decomposition can often be too restrictive, and so instead it is better to decompose `S` *approximately* as a low-rank matrix `T` plus a diagonal matrix `P`. 

### Our approach

This page contains several sample implementations of an estimation procedure for FA as described in [Bertsimas, Copenhaver, and Mazumder, "Certifiably Optimal Low Rank Factor Analysis", Journal of Machine Learning Research 18 (2017) ("BCM17")](http://jmlr.org/papers/v18/15-613.html). The approach is based on conditional gradient methods from convex optimization. We provide several sample implementations of Algorithm 1 (see page 13). Given the well-structured nature of the problems solved in Algorithm 1, there are algorithmic improvements that can be made to the implementations here, but these serve as a good starting point.

The problem that we focus on is the case of *q=1* outlined in BCM17. This special case is known as *Approximate Minimum Rank Factor Analysis*, or MRFA, and takes the form
```
minimize	nuclear_norm( S - ( T + P ) )
subject to 	P, T >= 0
		S - P >= 0
		P diagonal
		rank(T) <= r
```
Here the optimization variables are `T` and `P`; the notation `A >= 0` denotes that a matrix `A` is symmetric positive semidefinite; and the nuclear norm is the sum of the singular values. As shown in BCM17, this can be rewritten exactly as
```
minimize	trace( W*S - W*P )
subject to 	W, P >= 0
		S - P >= 0
		P diagonal
		I - W >= 0
		trace(W) = p - r
```
where `p` is the number of variables (i.e. `S` is a `p` by `p` matrix) and `I` is the `p` by `p` identity matrix.


The resulting algorithm described in BCM17 amounts to an alternating minimization scheme in `W` and `P`.


### Implementations

The two implementations are as follows:

1. [R](./R/)

   This implementation is as described in the paper. In particular, the update with respect to `P`, where `W` is fixed, is solved using the customized ADMM (Section 3.2.2).


2. [MATLAB](./matlab/)

   This implementation relies on [`cvx`](https://cvxr.com/cvx/ "CVX") and highlights the simplicity of the alternating minimization approach.



### Citation

If you would like to cite this work, please use the following citation for BCM17:
```
@article{bcm17,
  author  = {Dimitris Bertsimas and Martin S. Copenhaver and Rahul Mazumder},
  title   = {Certifiably Optimal Low Rank Factor Analysis},
  journal = {Journal of Machine Learning Research},
  year    = {2017},
  volume  = {18},
  number  = {29},
  pages   = {1-53},
  url     = {http://jmlr.org/papers/v18/15-613.html}
}
```

