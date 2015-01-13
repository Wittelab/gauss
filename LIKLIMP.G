/* needed to test: _bprior=0; _lp=0; _binit=0; */
/* This proc searches for penalized LR limits.
   1. CURRENTLY REQUIRES t2 FIXED AND DIAGONAL
      (algorithm not valid for nondiagonal t2)
   2. Not valid for those regressors penalized by a
      constraint b = zp*bs, where bs is unknown.
 Inputs: x = design matrix,
         y = vector of case counts,
         n = vector of totals (if scalar, all totals will be set to this n),
         is = index of coefficients to be modelled at second stage
              -- if const = 1, the intercept will NOT be counted in the
                 indices and will NOT be modelled
              -- if is = 0, all coefficients except the intercept
                 will be modelled
         zp = second-stage design matrix for linear constraint
              (set to 0 for shrinkage of coefficients to bp)
         t2 = second-stage residual variances
              -- set to 0 for single estimated tau2 (empirical Bayes),
         rep = vector of repetion counts (0 if none),
         const = 0 if no intercept (constant) term in model
                 (limits won't be computed for constant),
         offset = vector of offsets (0 if no offsets),
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
         yname = scalar name of outcome (may be 0 if no printed output).
         bnames2 = vector of second-stage coefficient names
                  (set to 0 if no second-stage printed output is desired),
         indl = vector of indices of x cols for which limits will be found
                (all if 0),
         bp = point to shrink to if zp is 0 (default is 0)
Outputs: b = coefficient estimates,
         bcov = inverse-information covariance matrix,
         bs = second-stage coefficient estimates (0 if zp eq 0),
         vs = covariance matrix for bs (0 if zp eq 0),
         t2 = t2 if input t2 gt 0, estimate of t2 if input t2 eq 0
         dev = residual deviance for model,
         devp = residual penalized deviance,
         df = residual degrees of freedom.
         limits = k by 2 matrix of lower and upper limits for b[indl].
Globals (may be set by user from calling program; otherwise defaults are used):
      All the globals for logregp.g, plus the following declared here:
      _clevel = confidence level expressed as proportion (default is .95)
      _criter = convergence criterion for limits
      _laxit = maximum number of iterations for limits (default is 30)
*/
proc (9)
     = liklimp(x,y,n,is,zp,t2,rep,const,offset,bnames,yname,bnames2,indl,bp);
declare matrix _criter = .002; declare matrix _laxit = 30;
declare matrix _clevel = .95;
local t0,np,snp,nis,snis,t2i,bhat,bcov,bs,vs,t2ret,dev,devp,dfp,se,b1,z,limits,
      j,ind,indis,irs,ips,iss,ii,jpen,b,lim,i,cnv,rad,offs,devpl,df,dp,radold;
t0 = date; bnames = 0$+bnames;
np = cols(x); snp = seqa(1,1,np);
if is[1] eq 0; is = snp; endif;
nis = rows(is); snis = seqa(1,1,rows(is));
if indl[1] eq 0; indl = snp; endif;
if cols(t2) eq 1; t2 = t2.*eye(nis); endif; t2i = invpd(t2);
zp = zp.*ones(nis,1); bp = bp.*ones(nis,1); _bprior = bp; _lp = 0;
if bnames[1]; print;
   "PROC LIKLIMP.G: penalized likelihood limits for logistic regression --";;
   print;
endif;
{ bhat,bcov,bs,vs,t2ret,dev,devp,dfp } =
  logregp(x,y,n,is,zp,t2,rep,const,offset,bnames,yname,bnames2);
se = sqrt(diag(bcov));
if const; b1 = bhat[1]; bhat = bhat[2:np+1]; se = se[2:np+1]; endif;
/* Desired Z value: */ z = cinvnorm((1-_clevel)/2);
/* Storage for limits: */ limits = zeros(rows(indl),2);
/* Loop over indl: */ j = 0;
do until j ge rows(indl); j = j + 1;
   /* Indices in b and in is of variable to be offset (for which limits are
      to be tested): */ ind = indl[j]; indis = snis'(snis .eq ind);
   /* Define vectors of regressor indices (excluding ind): */
      /* in bhat: */ irs = selif(snp,snp .ne ind);
      /* in bp & t2: */ ips = selif(snis,snis .ne indis);
      /* in is: */ iss = selif(seqa(1,1,np-1),sumr(irs .eq is'));
   /* resorted indices: */ ii = ind|irs;
   if bnames[1] /* print headers for this variable: */;
      format /rdn 3,0; print; 100*(1-cdfnt(z));;
      format 10,5; "% likelihood limits for "$+bnames[ind]$+":";
   endif;
   jpen = rows(iss) lt rows(is);
   b = bhat[irs]; if const; b = b1|_binit; endif;
   _bprior = bp[ips]; _lp = 0;
   /* Initialize CI radius at Wald radius: */ rad = z*se[ind];
   if bnames[1]; "   Wald lower limit for OR :";; exp(bhat[ind] - rad); endif;
   /* Iterate to lower limit using offset to find where deviance test yields
      desired z value: */
   i = 0; cnv = 0;
   do until (i gt _laxit) or cnv; i = i + 1;
      lim = bhat[ind] - rad; offs = lim*x[.,ind] + offset; _binit = b;
      if jpen; _lp = (lim - bp[indis])*t2i[ips,indis]; endif;
      { b,bcov,bs,vs,t2ret,dev,devpl,dfp } =
        logregp(x[.,irs],y,n,iss,zp[ips,.],t2[ips,ips],rep,const,offs,0,0,0);
      if jpen; dp = (lim|b[const+iss]) - bp[ii];
               devpl = dev + dp't2i[ii,ii]*dp; endif;
      radold = rad;
      rad = rad*z/sqrt(devpl-devp);
      /* Test for convergence: */ cnv = abs(rad - radold) lt _criter;
      if bnames[1] /* print results for this variable: */;
         format /rdn 3,0; "     At iteration ";; i;; format 10,5;
         ", lower limit & p = ";; exp(lim);; cdfchic(devpl-devp,1);
      endif;
   endo;
   if bnames[1]; "   Final lower limit for OR:";; exp(lim); print; endif;
   limits[j,1] = lim; rad = z*se[ind];
   if bnames[1]; "   Wald upper limit for OR :";; exp(bhat[ind] + rad); endif;
   b = bhat[irs]; if const; b = b1|_binit; endif;
   /* Iterate to upper limit: */
   i = 0; cnv = 0;
   do until (i gt _laxit) or cnv; i = i + 1;
      lim = bhat[ind] + rad; offs = lim*x[.,ind] + offset; _binit = b;
      if jpen; _lp = (lim - bp[indis])*t2i[ips,indis]; endif;
      { b,bcov,bs,vs,t2ret,dev,devpl,dfp } =
        logregp(x[.,irs],y,n,iss,zp[ips,.],t2[ips,ips],rep,const,offs,0,0,0);
      if jpen; dp = (lim|b[const+iss]) - bp[ii];
               devpl = dev + dp't2i[ii,ii]*dp; endif;
      radold = rad;
      rad = rad*z/sqrt(devpl-devp);
      cnv = abs(rad - radold) lt _criter;
      if bnames[1];
         format /rdn 3,0; "     At iteration ";; i;; format 10,5;
         ", upper limit & p = ";; exp(lim);; cdfchic(devpl-devp,1);
      endif;
   endo;
   if bnames[1]; "   Final upper limit for OR:";; exp(lim); endif;
   limits[j,2] = lim;
endo;
_bprior = bp; _lp = 0; _binit = bhat;
"Total run time: ";; etstr(ethsec(t0,date));
retp(b,bcov,bs,vs,t2,dev,devp,dfp,limits);
endp;
