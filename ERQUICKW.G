/* This proc computes ML Poisson-exponential regression coefficients with external weights u
   - no outputs or special options - returns only b, bcov, cnv.
 Inputs: x = design matrix (MUST INCLUDE CONSTANT IF ONE IS DESIRED),
         y = vector of case counts,
         n = vector of person-times,
         u = external weight,
         offset = vector of offsets (0 if no offsets),
         b = initial values for coefficients.
 Outputs: b = coefficient estimates (0 if singularity occurs),
          bcov = see sandwich covariance matrix option
          cnv = 1 if convergence in _maxit iterations, 0 if not,
                -1 if singularity occurs.
 Globals: _specrob = 0 for sandwich covariance assuming Poisson,
                                    1 for inverse Poisson information,
                                    2 for sandwich covariance assuming binomial.
*/
proc (3) = erquickw(x,y,n,u,offset,b);
/* Max. iterations    : */ declare matrix _maxit = 30;
/* Convergence criter : */ declare matrix _bcriter = .001;
/* Covariance estimate: */ declare matrix _specrob = 0;
local bold,bcov,iter,cnv,p,mu,ss;
/* Initialize beta, cnv, counter: */
   b = b.*ones(cols(x),1); cnv = 0; iter = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* save old values: */ bold = b;
   /* expected counts: */ mu = n.*exp(offset+x*b);
   /* covariance of b: */ trap 1; bcov = invpd(x'((u.*mu).*x));
   if scalerr(bcov); /* Warn & return if singularity occurs: */
      "SINGULARITY IN ERQUICKW.G AT ITERATION";; iter; retp(0,0,-1); endif; trap 0;
   /* Newton step: */ b = b + bcov*((x.*u)'(y-mu));
   cnv = prodc(abs(b - bold) .lt _bcriter);
endo;
if cnv le 0; "NONCONVERGENCE IN LRQUICKW"; endif; 
if _specrob eq 0; 
         ss = x.*u.*(y-mu);  bcov = bcov'(ss'ss)*bcov;
   elseif _specrob eq 2;
   bcov = bcov*(x'((u.*u.*(y + mu.*(mu - 2*y)./(n + n .eq 0))).*x))*bcov; 
endif;
retp(b,bcov,cnv);
endp;

