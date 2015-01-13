/* This proc does penalized-likelihood semi-Bayes exponential regression
   - no outputs or special options - returns only b,bcov,bs,vs,cnv.
 Inputs: x = design matrix (MUST INCLUDE CONSTANT IF ONE IS DESIRED),
         y = vector of case counts,
         n = vector of totals,
         z = second-stage design matrix for linear constraint
             (set to 0 for shrinkage of coefficients to _bp)
         t2 = second-stage residual variances
         offset = vector of offsets (0 if no offsets),
         b = initial values for coefficients,
 Outputs: b = coefficient estimates (0 if singularity occurs),
          bcov = inverse-information covariance matrix,
          bs = second-stage coefficient estimates (0 if z eq 0),
          vs = covariance matrix for bs (0 if z eq 0),
          cnv = 1 if convergence in _maxit iterations, 0 if not,
                -1 if singularity occurs.
  Globals: _bp = point to shrink to if z is 0
*/
proc (5) = erpquick(x,y,n,z,t2,offset,b);
/* Max. iterations: */ declare matrix _maxit = 30;
/* b convergence criterion: */ declare matrix _bcriter = .001;
/* t convergence criterion: */ declare matrix _tcriter = .01;
/* prior means if z = 0: */ declare matrix _bp = 0;
local np,eyenp,df2,tobp,ip,iter,cnv,bold,mu,
      inf,infi,w2,wz,ih,bcov,bs,vs,wvce;
np = cols(x); eyenp = eye(np);
tobp = sumall(abs(z)) eq 0;
if tobp; ip = eyenp./t2; endif;
   if rows(b) ne np;
      /* start with WLS estimates: */ local r;
      r = (y+1)./(n + sumc(n)/sumc(y)); mu = n.*r;
      trap 1; bcov = invpd(moment(sqrt(mu).*x,0));
      if scalerr(bcov); "SINGULARITY INITIALIZING ERPQUICK.G";
         retp(0,0,0,0,0,0,0,0,0); endif; trap 0;
      b = bcov*(x'(mu.*(ln(r)-offset)));
   endif;
/* Initialize: */ iter = 0; cnv = 0;
/* Begin iterative reweighted penalized LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* save old values: */ bold = b;
   /* expected counts: */ mu = n.*exp(offset+x*b);
   /* first-stage information matrix: */ inf = x'(mu.*x);
   trap 1; infi = invpd(inf);
      if scalerr(infi); /* Warn & return if singularity occurs: */
         "SINGULAR INFORMATION IN ERPQUICK.G "; retp(0,0,0,0,-1); endif;
      trap 0;
   if not tobp;
      /* inverse of estimated unconditional covariance of b, times z: */
         w2 = invpd(infi + t2.*eyenp); wz = w2*z;
      /* 2nd-stage I-H (residual projection) matrix */
         ih = eyenp - z*invpd(z'wz)*wz';
       /* prior information: */ ip = ih'(ih./t2);
   endif;
   /* inverse 2nd derivative of penalized LL: */
      trap 1; bcov = invpd(inf + ip); trap 0;
      if scalerr(bcov); /* Warn & return if singularity occurs: */
         "SINGULAR 2ND DERIVATIVE IN ERPQUICK.G"; retp(0,0,0,0,-1); endif;
   /* Newton step: */ b = b + bcov*(x'(y-mu) - ip*(b-_bp));
   cnv = prodc(abs(b - bold) .lt _bcriter);
endo;
if not cnv; "NO CONVERGENCE IN ERPQUICK.G"; endif;
/* second-stage covariance matrix, coefficients, and expected b: */
   if tobp; vs = 0; bs = 0; else; vs = invpd(z'wz); bs = vs*(wz'b); endif;
   if not tobp; bcov = infi - infi*w2*ih*infi;
      else; bcov = infi - infi*invpd(infi + t2.*eyenp)*infi;
   endif;
retp(b,bcov,bs,vs,cnv);
endp;

