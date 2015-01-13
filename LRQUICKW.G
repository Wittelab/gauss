/* This proc computes ML binomial-logistic regression coefficients with external weights u
   - no outputs or special options - returns only b, bcov, cnv.
 Inputs: x = design matrix (MUST INCLUDE CONSTANT IF ONE IS DESIRED),
         y = vector of case counts,
         n = vector of totals,
         u = external weight,
         offset = vector of offsets (0 if no offsets),
         b = initial values for coefficients.
 Outputs: b = coefficient estimates (0 if singularity occurs),
          bcov = see sandwich covariance matrix option
          cnv = 1 if convergence in _maxit iterations, 0 if not,
                -1 if singularity occurs.
 Globals: _specrob = 1 for sandwich covariance (default is 1), 0 for inverse information
*/
proc (3) = lrquickw(x,y,n,u,offset,b);
/* Max. iterations    : */ declare matrix _maxit = 30;
/* Convergence criter : */ declare matrix _bcriter = .001;
/* Covariance estimate: */ declare matrix _specrob = 1;
local bold,bcov,iter,cnv,p,mu;
/* Initialize beta, cnv, counter: */
   b = b.*ones(cols(x),1); cnv = 0; iter = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* save old values: */ bold = b;
   /* expected proportions: */ p = 1/(1 + exp(-offset-x*b));
   /* expected counts: */ mu = n.*p;
   /* covariance of b: */ trap 1; bcov = invpd(x'((u.*mu.*(1-p)).*x));
   if scalerr(bcov); /* Warn & return if singularity occurs: */
      "SINGULARITY IN LRQUICKW.G AT ITERATION";; iter; retp(0,0,-1); endif; trap 0;
   /* Newton step: */ b = b + bcov*((x.*u)'(y-mu));
   cnv = prodc(abs(b - bold) .lt _bcriter);
endo;
if cnv le 0; "NONCONVERGENCE IN LRQUICKW"; endif;
if _specrob; 
   bcov = bcov*(x'((u.*u.*(y.*(1-2*p) + mu.*p)).*x))*bcov; 
endif;
retp(b,bcov,cnv);
endp;

