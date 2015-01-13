/* This procedure computes the Breslow-Cain pseudo-maximum
   likelihood estimates for logistic regression analysis of two-stage
   case-control data. Note that a constant term is required and
   automatically added to the model. See Breslow, N.E. and Cain, K.C.
   (1988), Biometrika, v. 75, p. 11-20.
 Inputs:
   1) SECOND-STAGE SAMPLE DATA:
        x = regressor matrix (may include effects of
            first-stage variables used to construct st).
        y = vector of case counts.
        n = vector of totals (if scalar, all totals will be set to this n).
        st = discrete vector giving the second-stage sampling
             stratum to which each observation belongs.
   2) FIRST-STAGE SAMPLE DATA:
        st1 = discrete vector listing the second-stage sampling strata
              (this lists the unique values in st).
        y1 = vector of first-stage case counts for st1.
        z1 = vector of first-stage control (noncase) counts for st1.
   3) Miscellaneous inputs:
        fixed = 1 if 2nd-stage sampling fractions are fixed, 0 if random
        rep = vector of repetion counts for 2nd-stage sample (0 if none).
        bnames = vector of coefficient names
                 (set to 0 if no printed output is desired).
        yname = scalar name of outcome (may be 0 if no printed output).
 Outputs: b = coefficient estimates,
          bcov = Breslow & Cain's preferred covariance matrix estimate
                 (may be nonpositive),
          bcovr = robust ("sandwich") covariance matrix for b,
          dev = residual pseudo-deviance for model,
          rss = residual weighted sum of squares,
          df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
*/
proc (6) = lr2stage(x,y,n,st,st1,y1,z1,fixed,rep,bnames,yname);
declare matrix _binit = 0; declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .1;
local t0,nr,np,y2,z2,offset,llsat,b,bold,infi,dev,devold,df,neq0,iter,cnv,
      eta,p,mu,w,xws,aba,bcov,s,ssy,ssz,g,bcovr,rss,mcc;
t0 = date;
if bnames[1]; print;
   " Proc lr2stage.g: logistic regression for two-stage case-control data";
   "      using the method of Breslow & Cain (Biometrika, 1988)."; print;
   if rows(bnames) ne cols(x);
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
endif;
nr = rows(x);
if rows(n) .eq 1; n = n*ones(nr,1); endif;
if sumc(rep) ne 0; y = y.*rep; n = n.*rep; endif;
/* Include constant (required: */  x = ones(nr,1)~x;
/* No. parameters: */ np = cols(x);
/* No. empty covariate patterns: */ neq0 = n .eq 0;
/* Degrees of freedom: */ df = nr - np - sumc(neq0);
/* Indicators of stratum membership: */ st = st .eq st1';
/* Sizes of 2nd-stage case and control samples: */ y2 = st'y; z2 = st'n - y2;
if sumc((y2|z2) .le 0); "ERROR";
 "There must be at least one case and one control in each 2nd-stage stratum.";
 "Case & control counts by stratum:";; y2~z2; end;
endif;
if sumc((y2|z2) .gt (y1|z1)); "ERROR";
 "There can be no more subjects in 2nd-stage than in 1st stage.";
 "Case & control counts by stage & stratum:";; y1~z1~y2~z2; end;
endif;
/* Offsets: */ offset = st*ln(y2.*z1./(z2.*y1));
/* Saturated pseudologlikelihood + n'ln(n) - sumc(ln(combin(n,y))): */
   llsat = y'ln(y + (y.==0)) + (n - y)'ln(n - y + (y.==n));
/* Initialize beta, deviance, convergence indicator, counter: */
   if rows(_binit) eq np;
      b = _binit;
      else; mu = y + sumc(y)/sumc(n); p = mu./(n+1);
      w = n.*p.*(1-p); b = invpd(x'(w.*x))*((w.*x)'(logit(p)-offset));
   endif;
   dev = 0; cnv = 0; iter = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictors: */ eta = offset + x*b;
   /* expected proportions: */ p = 1./(1 + exp(-eta));
   /* expected counts: */ mu = n.*p;
   devold = dev; dev = 2*(llsat - (y'eta + n'ln(n + neq0 - mu)));
   if bnames[1]; "Pseudo-deviance at iteration ";;
      format /rds 4,0; iter;; ":";; format 12,3; dev; endif;
   if iter gt 1 and devold lt dev; dev = devold;
      /* step halve: */ b = (bold + b)/2; continue;
   endif;
   /* weight vector: */ w = mu.*(1-p);
   /* inverse information: */ infi = invpd((w.*x)'x);
   /* Newton step: */ bold = b; b = b + infi*(x'(y-mu));
   cnv = prodc(abs(b - bold) .lt _bcriter) and (abs(devold - dev) lt _dcriter);
endo;
if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
/* Breslow-Cain "W" vectors: */ xws = (x.*w)'st;
   if fixed; aba = 0; else;
      /* include est. of variability from sampling-fraction variation: */
         aba = xws*((1/y1 + 1/z1).*xws')
               - sumr(xws)*(1/sumc(y1) + 1/sumc(z1))*sumr(xws)';
   endif;
/* Breslow-Cain covariance est: */
   bcov = infi - infi*(xws*((1/y2 + 1/z2).*xws') - aba)*infi;
/* score contributions: */ s = x.*(y - mu);
/* case & control statum-specific score sums: */
   ssy = ((y .eq n).*s)'st; ssz = ((y .eq 0).*s)'st;
/* Breslow-Cain G-matrix est.: */
   g = (s'(s./(n + neq0)) - (ssy./y2')*ssy' - (ssz./z2')*ssz');
/* "robust" covariance est.: */ bcovr = infi*(g + aba)*infi;
/* weighted sum of squared residuals: */ rss = (y-mu)'((y-mu)./(w + neq0));
if bnames[1]; print; mcc = meanc(y)|meanc(n-y);
   if fixed; "Variances computed using fixed sampling-fraction formula.";
      else; "Variances computed using random sampling fraction formula.";
   endif; print;
   if prodc(diag(bcov) gt 0); "With Breslow-Cain (1988, eq. 9) variances:";;
      rreport(b,sqrt(diag(bcov)),bnames,yname,-1,-1,0,1); print;
      else; "Breslow-Cain (1988, eq. 9) variances not all positive.";
   endif;
   "Using specification-robust variances:";;
   rreport(b,sqrt(diag(bcovr)),bnames,yname,dev,rss,df,1+(minc(mcc) ge 3));
   "The means of the 2nd-stage case & noncase counts are ";; format 6,1; mcc';
   if minc(mcc) lt 5;
         "The minimum of the mean case & mean noncase counts is less than 5.";
      if minc(mcc) ge 3;
         "Hence the above tests of fit may be invalid.";
      endif;
      if sumc(n gt 2) gt nr/2 and df gt 50;
         " Most of the denominators are 2 or more and df > 50. For this";
         " situation, the following conditional approximate Gaussian test";
         " of fit is suggested by McCullagh and Nelder (1989, p.122):";
         local v,wv,q,xq,cxq,erss,vrss,zrss;
         v = sumr((x*bcov).*x);
         wv = w.*v;
         q = 1 - 2*mu./n + neq0;
         xq = (q'x)';
         cxq = bcov*xq;
         erss = df + 3*sumc(wv./n + neq0) + (((wv.*q)'x)*cxq - sumc(v))/2;
         vrss = (1 - np/nr)*(2*(nr - sumc(1/n + neq0)) + xq'cxq);
         zrss = (rss - erss)/sqrt(vrss);
         " Approximate condtional E(RSS), Var(RSS),";;
         " Z-statistic, and upper p-value:"; format 12,3;
         "                 ";; erss;; vrss;; zrss;; cdfnc(zrss);
      else;
        "Because of the small counts in these data, one must resort to a";
        " grouped statistic (e.g., Hosmer-Lemeshow) or an embedding-model";
        " statistic (e.g., Brown's score statistic) to test fit.";
      endif;
   endif;
   "Total run time: ";; etstr(ethsec(t0,date));
endif; format /rds 8,3;
retp(b,bcov,bcovr,dev,rss,df);
endp;
