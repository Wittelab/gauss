/* This proc computes ML exponential-Poisson regression coefficients
   - no outputs or special options - returns only b, bcov, cnv.
 Inputs: x = design matrix (MUST INCLUDE CONSTANT IF ONE IS DESIRED),
         y = vector of case counts,
         n = vector of person-time (1 if none),
         offset = vector of offsets (0 if no offsets),
         b = initial values for coefficients.
 Outputs: b = coefficient estimates,
          bcov = inverse-information covariance matrix,
          iter = no. of iterations (no convergence if iter is negative).
*/
proc (3) = erquick(x,y,n,offset,b);
/* Max. iterations: */ declare matrix _maxit = 30;
/* Convergence criterion: */ declare matrix _bcriter = .001;
local bold,bcov,iter,cnv,mu;
/* Initialize beta, cnv, counter: */
   b = b.*ones(cols(x),1); cnv = 0; iter = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* save old values: */ bold = b;
   /* expected counts: */ mu = n.*exp(offset + x*b);
   trap 1; bcov = invpd(moment(sqrt(mu).*x,0));
      if scalerr(bcov); "SINGULARITY IN ERQUICK.G AT ITERATION";; iter;
         retp(b,0,-iter); endif; trap 0;
   /* Newton step: */ b = b + bcov*(x'(y-mu));
   cnv = prodc(abs(b - bold) .lt _bcriter);
endo;
if iter ge _maxit; 
    "NO CONVERGENCE IN ERQUICK.G AT ITERATION";; iter; iter = -iter; endif;
retp(b,bcov,iter);
endp;
