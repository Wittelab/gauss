/* This proc does penalized-likelihood logistic regression
   - no outputs or special options - returns only b,bcov,bs,vs,t2,cnv.
 Inputs: x = design matrix (MUST INCLUDE CONSTANT IF ONE IS DESIRED),
         y = vector of case counts,
         n = vector of totals,
         z = second-stage design matrix for linear constraint
             (set to 0 for shrinkage of coefficients to _bp)
         t2 = second-stage residual variances
              -- set to 0 for single estimated tau2 (empirical Bayes),
         offset = vector of offsets (0 if no offsets),
         b = initial values for coefficients,
         The following input values are ignored unless EB estimation is used:
             t2init = initial value for t2
             bc = 1 if t2 to be estimated by Breslow-Clayton (1993) formula
                  (which runs very slowly and takes much memory but is more
                  accurate in small samples), 0 for Morris estimate.
 Outputs: b = coefficient estimates (0 if singularity occurs),
          bcov = inverse-information covariance matrix,
          bs = second-stage coefficient estimates (0 if z eq 0),
          vs = covariance matrix for bs (0 if z eq 0),
          t2 = t2 if input t2 gt 0, estimate of t2 if input t2 eq 0
          cnv = 1 if convergence in _maxit iterations, 0 if not,
                -1 if singularity occurs.
  Globals: _bp = point to shrink to if z is 0
*/
proc (6) = lrpquick(x,y,n,z,t2,offset,b,t2init,bc);
/* Max. iterations: */ declare matrix _maxit = 30;
/* b convergence criterion: */ declare matrix _bcriter = .001;
/* t convergence criterion: */ declare matrix _tcriter = .01;
/* prior means if z = 0: */ declare matrix _bp = 0;
local np,eyenp,eb,df2,tobp,ip,iter,cnv,undsp,bold,t2old,p,mu,w,
      inf,infi,w2,wz,ih,e2,sw2,bcov,bs,vs,wvce;
np = cols(x); eyenp = eye(np);
eb = sumall(t2) eq 0;
if eb; t2 = t2init; df2 = np - cols(z);
   if bc; /* for Breslow-Clayton variance: */
      local xx,dxx,u,vi,xvu,vix; xx = x*x'; dxx = diag(xx); endif;
endif;
tobp = sumall(abs(z)) eq 0;
if tobp; ip = eyenp./t2; endif;
/* Initialize: */ b = b.*ones(cols(x),1); cnv = 0; undsp = 0; t2old = t2; iter = 0;
/* Begin iterative reweighted penalized LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* save old values: */ bold = b;
   /* expected proportions: */ p = 1/(1 + exp(-offset-x*b));
   /* expected counts: */ mu = n.*p;
   /* 1st-stage weight: */ w = mu.*(1-p);
   /* first-stage information matrix: */ inf = x'(w.*x);
   trap 1; infi = invpd(inf);
      if scalerr(infi); /* Warn & return if singularity occurs: */
         "SINGULAR INFORMATION IN LRPQUICK.G AT ITERATION";; iter; retp(0,0,0,0,t2,-1);
      endif;
      trap 0;
   if not tobp;
      /* inverse of estimated unconditional covariance of b, times z: */
         w2 = invpd(infi + t2.*eyenp); wz = w2*z;
      /* 2nd-stage I-H (residual projection) matrix */
         ih = eyenp - z*invpd(z'wz)*wz';
      if eb; /* update t2: */ t2old = t2;
         if bc; /* update Breslow-Clayton variance estimate: */
             u = (x*(ih*b) + (y - mu)./w);
             vi = invpd(diagrv(xx,dxx + 1/(w*t2)));
             xvu = x'(vi*u); vix = vi*xx;
             t2 = t2 + (np*(xvu'xvu)/df2 - t2*tracemat(vix))/sumall(vix.*vix');
             else; /* update Morris variance estimate: */
                   /* 1-step approx. to MLE: bml1 = b + infi*(x'(y-mu)); */
                   /* 1-step approx to 2nd-stage res, bml - (eyenp-ih)*b: */
                      e2 = infi*(x'(y-mu)) + ih*b; sw2 = sumall(w2);
                   t2 = (np*(e2'(w2*e2))/df2 - sumall(w2*infi))/sw2;
          endif;
          if t2 le .0001; "UNDERDISPERSION IN LRPQUICK.G";
             eb = 0; undsp = 1; t2 = .0001; t2old = t2; endif;
       endif;
       /* prior information: */ ip = ih'(ih./t2);
   endif;
   /* inverse 2nd derivative of penalized LL: */
      trap 1; bcov = invpd(inf + ip); trap 0;
      if scalerr(bcov); /* Warn & return if singularity occurs: */
         "SINGULAR 2ND DERIVATIVE IN LRPQUICK.G"; retp(0,0,0,0,t2,-1); endif;
   /* Newton step: */ b = b + bcov*(x'(y-mu) - ip*(b-_bp));
   cnv = prodc(abs(b - bold) .lt _bcriter)
         and (abs(sqrt(t2) - sqrt(t2old)) lt _tcriter);
endo;
if not cnv; "NO CONVERGENCE IN LRPQUICK.G"; endif;
/* second-stage covariance matrix, coefficients, and expected b: */
   if tobp; vs = 0; bs = 0; else; vs = invpd(z'wz); bs = vs*(wz'b); endif;
if eb or undsp; bcov = infi - ((df2-2)/df2)*infi*w2*ih*infi;
   if bc; e2 = infi*(x'(y-mu)) + ih*b; sw2 = sumall(w2); endif;
   /* add variance component to account for estimation of t2: */
   wvce = ((df2-2)/df2)*(w2*(infi*e2));
   wvce = wvce.*sqrt((sumall(w2*infi)/sw2 + t2)./(diag(infi) + t2));
   bcov = bcov + (2/(df2-2))*wvce*wvce';
   elseif not tobp; /* fixed-tau2 cov: */ bcov = infi - infi*w2*ih*infi;
   else; /* fixed-tau2 cov: */ bcov = infi - infi*invpd(infi + t2.*eyenp)*infi;
endif;
retp(b,bcov,bs,vs,t2,cnv);
endp;

