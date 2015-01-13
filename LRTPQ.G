/* This proc computes truncated-Poisson logistic regression coefficients
   - no outputs or special options - returns only b, bcov, cnv.
 Inputs: x = design matrix (MUST INCLUDE CONSTANT IF ONE IS DESIRED),
         y = vector of truncated Poisson counts
         m = maximum value for y (scalar)
         rep = repetition count (set to 1 for no count)
         offset = vector of offsets (0 if no offsets),
         binit = initial values for coefficients.
 Outputs: b = coefficient estimates,
          bcov = inverse-information covariance matrix,
          cnv = 1 if convergence in _maxit iterations, 0 if not.
*/
proc (3) = lrtpq(x,y,m,rep,offset,binit);
/* Max. iterations: */ declare matrix _maxit = 30;
/* Convergence criterion: */ declare matrix _bcriter = .001;
local m2,b,bold,bcov,iter,cnv,p,mu;
m = seqa(1,1,m); m2 = m.*m;
offset = offset - lnfact(m)';
/* Initialize beta, loglikelihood, counter: */
   b = binit; iter = 0; cnv = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit);
   iter = iter + 1;
   /* save old values: */ bold = b;
   p = exp(offset + (x*b).*m');
   /* expected proportions: */ p = p./(1 + sumr(p));
   /* expected counts: */ mu = p*m;
   /* covariance of b: */ bcov = invpd(x'((rep.*(p*m2 - mu.*mu)).*x));
   /* Newton step: */ b = b + bcov*(x'(rep.*(y-mu)));
   cnv = prodc(abs(b - bold) .lt _bcriter);
endo;
retp(b,bcov,cnv);
endp;
