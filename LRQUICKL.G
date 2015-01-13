/* This proc computes ML binomial-logistic regression coefficients
   - no outputs or special options - returns only b, bcov, llik, cnv.
 Inputs: x = design matrix (MUST INCLUDE CONSTANT IF ONE IS DESIRED),
         y = vector of case counts,
         n = vector of totals,
         offset = vector of offsets (0 if no offsets),
         b = initial values for coefficients.
 Outputs: b = coefficient estimates (0 if singularity occurs),
          bcov = inverse-information covariance matrix,
          llik = loglikelihood at convergence,
          cnv = 1 if convergence in _maxit iterations, 0 if not,
                -1 if singularity or divergence occurs.
*/
proc (4) = lrquickl(x,y,n,offset,b);
/* Max. iterations: */ declare matrix _maxit = 30;
/* Convergence criterion: */ declare matrix _bcriter = .001;
local bold,bcov,iter,cnv,eta,p,mu,sflag,llik;
/* Initialize beta, cnv, counter: */
   b = b.*ones(cols(x),1); bold = b; cnv = 0; iter = 0;
/* Flag for step fractioning: */ sflag = 4;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictor: */ eta = offset + x*b;
   /* expected proportions: */ p = 1/(1 + exp(-eta));
   /* expected counts: */ mu = n.*p;
   bcov = x'((mu.*(1-p)).*x);
   /* covariance of b: */ trap 1; bcov = invpd(bcov);
   if scalerr(bcov);
      if sflag; sflag = sflag-1;
         if iter eq 1; bold = 0; endif;
         /* reduce Newton step: */ b = .8*bold+.2*b; trap 0; continue;
         else; /* Warn & return if singularity occurs repeatedly: */
               "SINGULARITY IN LRQUICKL.G AT ITERATION";; iter; trap 0; retp(0,0,0,-1);
      endif;
   endif; trap 0;
   /* save old values: */ bold = b;
   /* Newton step: */ b = b + bcov*(x'(y-mu));
   cnv = prodc(abs(b - bold) .lt _bcriter);
endo;
if not cnv; "NONCONVERGENCE IN LRQUICKL.G"; endif;
/* loglikelihood + constant: */ llik = y'eta + n'ln(n-mu+((n-mu) .eq 0));
retp(b,bcov,llik,cnv);
endp;

