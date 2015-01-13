/* This procedure does ML logistic regression
 Inputs: x = design matrix,
         y = vector of case counts,
         n = vector of totals (if scalar, all totals will be set to this n),
         rep = vector of repetion counts (0 if none),
         const = 1 if you wish the program to add a constant to x, 0 if not;
                 if the first j (>1) columns of x represent a complete set
                 of indicators for one covariate, const = j will result in
                 weighted-mean centering of the remaining cols(x)-j covariates,
                 where the weights are those obtained from the final iteration,
                 which "approximately diagonalizes" the covariance submatrix
                 of the j indicator-coefficient (intercept) estimates.
         offset = vector of offsets (0 if no offsets),
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
         yname = scalar name of outcome (may be 0 if no printed output).
 Outputs: b = coefficient estimates (with intercept first if const = 1)
          bcov = inverse-information covariance matrix estimate,
          bcovr = Huber robust ("sandwich") covariance matrix for b
                  -- this estimate does not depend on the binomial
                     variance assumption used to derive bcov
          dev = residual deviance for model,
          rss = residual weighted sum of squares,
          df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _collaps = 1 if data rows with identical covariates should be collapsed
                 together before fitting (default is 0)
      _clevel = confidence level expressed as a proportion (default is .95)
      _fit = 1 for iteration reporting and tests of fit, 0 if not (default is
             1, but no tests provided if mean of case or noncase count is <3)
      _mis = 0 if no missing values, 1 if complete-case analysis
      _mcode = missing-value scalar or column vector with element
               for each column of x (default -99999; ignored if _mis = 0)
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .01)
      _specrob = 1 if you want an additional version of output
                 in which the Huber specification-robust "sandwich"
                 covariance matrix estimate (bcovr) is used instead of
                 the usual inverse-information estimate (default is 0)
      _bias = 1 if you want an additional printout of odds ratios
              computed with the Cox-Snell ML bias correction (default is 0)
      _bprior = prior means (used for Bayesian estimation; default is 0)
*/
proc (6) = logreg(x,y,n,rep,const,offset,bnames,yname);
declare matrix _collaps = 0; declare matrix _clevel = .95;
declare matrix _mis = 0; declare matrix _mcode = -99999;
declare matrix _binit = 0; declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .01;
declare matrix _fit = 1; declare matrix _specrob = 0;
declare matrix _bias = 0; declare matrix _bprior = 0;
local t0,neq0,nr,mcc,np,llsat,b,bold,bcov,dev,devold,df,iter,cnv,
      eta,p,mu,w,rss,bcovr;
t0 = date; bnames = 0$+bnames; yname = 0$+yname;
if sumc(rep); y = y.*rep; n = n.*rep; else; n = n.*ones(rows(x),1); endif;
/* Delete empty covariate patterns: */ neq0 = n .eq 0;
   if _mis eq 1; /* also delete records with missing regressor values: */
      neq0 = neq0 .or sumr(x .eq _mcode'); endif;
   if sumc(neq0); x = delif(x,neq0); y = delif(y,neq0); n = delif(n,neq0);
      if rows(offset) eq rows(neq0); offset = delif(offset,neq0); endif;
   endif;
if _collaps; np = seqa(3,1,cols(x));
   if bnames[1]; "Collapsing over rows with identical covariates."; endif;
   if rows(offset) eq 1; x = matcol(y~n~x,np);
      else; x = matcol(y~n~x~offset,np); offset = x[.,cols(x)];
   endif;
   y = x[.,1]; n = x[.,2]; x = x[.,np];
endif;
nr = rows(x); np = cols(x);
if bnames[1]; print; "PROC LOGREG.G: ML logistic regression for "$+yname$+".";
   if rows(bnames) ne np;
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
   mcc = meanc(y)|meanc(n-y);
   print; "Nos. of cases and noncases :";; format /rds 6,0; rows(x)*mcc';
   "No. records used: ";; nr; "No. parameters  : ";; (const eq 1)+np;
endif;
if rows(_bprior) eq np; /*offset by prior*/ offset = offset + x*_bprior; endif;
if const eq 1; /* add constant: */ x = ones(nr,1)~x; np = 1+np; endif;
if rank(x'x) lt np; "DESIGN MATRIX RANK DEFICIENT IN LOGREG.G";
   retp(0,0,0,0,0,0); endif;
/* Degrees of freedom: */ df = nr - np;
/* Saturated loglikelihood + n'ln(n) - sumc(ln(combin(n,y))): */
   llsat = y'ln(y + (y.==0)) + (n - y)'ln(n - y + (y.==n));
/* Initialize beta, deviance, convergence indicator, counter: */
   if rows(_binit) eq np; b = _binit;
      else; /* start with WLS estimates: */
      p = (y + sumc(y)/sumc(n) + .01)./(n + 1.02);
      w = n.*p.*(1-p);
      trap 1; bcov = invpd(moment(sqrt(w).*x,0));
      if scalerr(bcov); "SINGULARITY INITIALIZING LOGREG.G";
         retp(0,0,0,0,0,0); endif; trap 0;
      b = bcov*(x'(w.*(logit(p)-offset)));
   endif;
dev = 0; cnv = 0; iter = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictors: */ eta = offset + x*b;
   /* expected proportions: */ p = 1./(1 + exp(-eta));
   if sumc(p .ge 1) or sumc(p .le 0); /* step halve: */
      b = (bold + b)/2; continue; endif;
   /* expected counts: */ mu = n.*p;
   devold = dev; dev = 2*(llsat - (y'eta + n'ln(n - mu)));
   if _fit*bnames[1]; "Deviance at iteration ";;
      format /rds 4,0; iter;; ":";; format 12,3; dev; endif;
   if iter gt 1 and devold lt dev; dev = devold;
      /* step halve: */ b = (bold + b)/2; continue;
   endif;
   /* weight vector: */ w = mu.*(1-p);
   /* covariance of b: */ trap 1; bcov = invpd(moment(sqrt(w).*x,0));
      if scalerr(bcov); "INFORMATION NOT POSITIVE DEFINITE IN LOGREG.G";
         retp(0,0,0,0,0,0); endif; trap 0;
   /* Newton step: */ bold = b; b = b + bcov*(x'(y-mu));
   cnv = prodc(abs(b - bold) .lt _bcriter) and (abs(devold - dev) lt _dcriter);
endo;
if iter ge _maxit; "WARNING: No convergence after ";;
   format 4,0; iter;; " iterations"; format 9,3;
endif;
/* add prior back in: */ b[1+const:np] = b[1+const:np] + _bprior;
/* weighted sum of squared residuals: */ rss = (y-mu)'((y-mu)./w);
/* pseudo-residual: e = (y - mu)./w; */
/* pseudo-observed logit: z = eta + e; */
/* note that b = bcov*wx'z; */
if const gt 1 and np gt const; local k; k = seqa(const+1,1,np-const);
   /* center all covariates after the first j, then recompute bcov & b: */
   x[.,k] = x[.,k] - (w/sumc(w))'x[.,k];
   bcov = invpd(moment(sqrt(w).*x,0)); b = bcov*(x'(w.*eta + y - mu));
endif;
/* the squared residuals are y.*(1-p)^2 + (n-y).*p^2 = y + p.*(mu - 2*y) */
/* hence the robust ("sandwich") covariance is: */
   bcovr = bcov*(x'((y.*(1-2*p)+ mu.*p).*x))*bcov; 
if bnames[1];
   if _specrob; "With inverse-information variances:"; endif;
   rreport(b,sqrt(diag(bcov)),bnames,yname,dev,
           rss*(-1)^(1-_fit),df,1+_fit*(minc(mcc) ge 3));
   if _fit;
      "The means of the case & noncase counts are ";; format 6,1; mcc';
      if minc(mcc) lt 5;
         "The minimum of the mean case & mean noncase counts is less than 5.";
         if minc(mcc) ge 3;
           "Hence the above deviance and rss tests of fit may be invalid.";
         endif;
         if sumc(n gt 2) gt nr/2 and df gt 50;
            " Most of the denominators are 2 or more and df > 50. For this";
            " situation, the following conditional approximate Gaussian test";
            " of fit is suggested by McCullagh and Nelder (1989, p.122):";
            local v,wv,q,xq,cxq,erss,vrss,zrss;
            v = sumr((x*bcov).*x);
            wv = w.*v;
            q = 1 - 2*mu./n;
            xq = (q'x)';
            cxq = bcov*xq;
            erss = df + 3*sumc(wv./n) + (((wv.*q)'x)*cxq - sumc(v))/2;
            vrss = (1 - np/nr)*(2*(nr - sumc(1/n)) + xq'cxq);
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
   endif;
   if _specrob; print; "With specification-robust variances:";
      rreport(b,sqrt(diag(bcovr)),bnames,yname,-1,-1,df,1); endif;
   if _bias; "Odds ratios:  ";;
       "ML, bias-corrected (& confidence limits), bias factors:";
       local bcovf,bfold,bf,rad,i;
       bcovf = bcov; bfold = 0; cnv = 0; iter = 0;
       bf = bcovf*(x'(w.*sumr((x*bcovf).*x).*(p - 1/2)));
       /* Begin iterative bias correction: */
       do until cnv or (iter ge _maxit); iter = iter + 1;
          eta = offset + x*(b - bf); p = 1./(1 + exp(-eta)); w = n.*p.*(1-p);
          trap 1; bcovf = invpd(moment(sqrt(w).*x,0));
          if scalerr(bcov); "SINGULARITY IN LOGREG.G BIAS CORRECTION";
             retp(0,0,0,0,0,0); endif; trap 0;
          bfold = bf; bf = bcovf*(x'(w.*sumr((x*bcovf).*x).*(p - 1/2)));
          cnv = prodc(abs(bf - bfold) .lt _bcriter);
       endo;
       if iter ge _maxit; "WARNING: No convergence of bias correction after ";;
          format 4,0; iter;; " iterations"; endif;
       format 10,2; rad = cdfni(1-(1-_clevel)/2).*((-1)~1); i = 0;
       do until i ge np; i = i + 1;
          $bnames[i];; exp(b[i]~(b[i]-bf[i]));;
          "(";;exp((b[i]-bf[i])+rad*sqrt(bcovf[i,i]));;")";; exp(bf[i]);
       endo;
       "Deviances  ";; dev;; 2*(llsat - (y'eta + n'ln(n - n.*p)));
   endif;
   if _fit; "Total run time: ";; etstr(ethsec(t0,date)); endif;
   print; format 10,3;
endif;
retp(b,bcov,bcovr,dev,rss,df);
endp;

