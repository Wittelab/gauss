/* This proc computes ML binomial-logistic regression coefficients
   - no outputs or special options - returns only b, bcov, cnv.
 Inputs: x = design matrix (MUST INCLUDE CONSTANT IF ONE IS DESIRED),
         y = vector of case counts,
         n = vector of totals,
         offset = vector of offsets (0 if no offsets),
         b = initial values for coefficients.
 Outputs: b = coefficient estimates (0 if singularity occurs),
          bcov = inverse-information covariance matrix,
          cnv = 1 if convergence in _maxit iterations, 0 if not,
                -1 if singularity occurs.
*/
proc (3) = lrquick(x,y,n,offset,b);
/* Max. iterations: */ declare matrix _maxit = 30;
/* Convergence criterion: */ declare matrix _bcriter = .001;
local bold,bcov,iter,cnv,p,mu;
/* Initialize beta, cnv, counter: */
   b = b.*ones(cols(x),1); cnv = 0; iter = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* save old values: */ bold = b;
   /* expected proportions: */ p = 1/(1 + exp(-offset-x*b));
   /* expected counts: */ mu = n.*p;
   /* covariance of b: */ trap 1; bcov = invpd(x'((mu.*(1-p)).*x));
   if scalerr(bcov); /* Warn & return if singularity occurs: */
      "SINGULARITY IN LRQUICK.G AT ITERATION";; iter; retp(b,0,-iter); endif; trap 0;
   /* Newton step: */ b = b + bcov*(x'(y-mu)); 
   cnv = prodc(abs(b - bold) .lt _bcriter); 
endo;
if iter ge _maxit; 
    "NO CONVERGENCE IN LRQUICK.G AT ITERATION";; iter; iter = -iter; endif;
retp(b,bcov,iter);
endp;

