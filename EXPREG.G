/* This procedure does ML exponential-Poisson regression
 Inputs: x = design matrix,
         y = vector of case counts,
         n = vector of person-time (if a scalar, all n will have this value)
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
                  (set to zero if no printed output is desired),
         yname = scalar name of outcome (may be 0 if no printed output).
 Outputs: b = coefficient estimates,
          bcov = inverse-information covariance matrix,
          bcovr = robust ("sandwich") covariance matrix for b
                  -- default assumes Poisson variation;
                  set _specrob = 2 for binomial data
          dev = residual deviance for model,
          rss = residual weighted sum of squares,
          df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _collaps = 1 if data rows with identical covariates should be collapsed
      _clevel = confidence level expressed as a proportion (default .95)
      _fit = 1 for iteration reporting and tests of fit, 0 if not (default is
             1, but no tests provided if mean case count is <3)
      _mis = 0 if no missing values, 1 if complete-case analysis
      _mcode = missing-value scalar or column vector with element
               for each column of x (default is -99999; ignored if _mis = 0)
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
      _specrob = 1 or 2 if you want an additional version of output
                 in which a structural-specification-robust "sandwich"
                 covariance matrix estimate (bcovr) is used instead of
                 the usual inverse-information estimate (default is 0):
                 1 assuming Poisson, 2 assuming binomial variability
                 (the binomial option should be used when expreg.g is used
                  for pseudo-likelihood fitting of exponential risk models).
*/
proc (6) = expreg(x,y,n,rep,const,offset,bnames,yname);
declare matrix _collaps = 0; declare matrix _clevel = .95;
declare matrix _mis = 0; declare matrix _mcode = -99999;
declare matrix _binit = 0; declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .1;
declare matrix _specrob = 0; declare matrix _fit = 1;
local t0,neq0,mcc,nr,np,llsat,b,r,bold,bcov,dev,devold,df,cnv,
      iter,eta,mu,rss,bcovr;
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
if bnames[1];
   print; " PROC EXPREG.G: ML Poisson exponential regression for "$+yname$+".";
   if rows(bnames) ne np;
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
   mcc = meanc(y); "Total no. of events (cases) and person-time: ";;
   format /rds 6,0; nr*mcc;; format 9,2; sumc(n); format 6,0;
   "No. records used: ";; nr; "No. parameters  : ";; (const eq 1)+np;
endif;
/* Add constant if requested: */
   if const eq 1; x = ones(nr,1)~x; np = 1 + np; endif;
if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN EXPREG.G";
   retp(0,0,0,0,0,0); endif;
/* Degrees of freedom: */ df = nr - np;
/* Saturated loglikelihood: */
   llsat = y'ln(y./n + (y.==0)) - sumc(y) /* - sumc(y!) */ ;
/* Initialize beta, deviance, convergence indicator, counter: */
   if rows(_binit) eq np; b = _binit;
      else; /* start with WLS estimates: */
      r = (y+1)./(n + sumc(n)/sumc(y)); mu = n.*r;
      trap 1; bcov = invpd(moment(sqrt(mu).*x,0));
      if scalerr(bcov); "SINGULARITY INITIALIZING EXPREG.G";
         retp(0,0,0,0,0,0); endif; trap 0;
      b = bcov*(x'(mu.*(ln(r)-offset)));
   endif;
   dev = 0; cnv = 0; iter = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictor: */ eta = offset + x*b;
   /* fitted case counts: */ mu =  n.*exp(eta);
   /* update deviance: */ devold = dev; dev = 2*(llsat - (y'eta - sumc(mu)));
   if _fit*bnames[1]; "Deviance at iteration ";;
      format /rds 4,0; iter;; ":";; format 12,3; dev; endif;
   if iter gt 1 and devold lt dev; dev = devold;
      /* step halve: */ b = (bold + b)/2; continue;
   endif;
   /* mu is weight vector in this case */
   /* covariance of b: */ trap 1; bcov = invpd(moment(sqrt(mu).*x,0));
      if scalerr(bcov); "INFORMATION NOT POSITIVE DEFINITE IN EXPREG.G";
         retp(0,0,0,0,0,0); endif; trap 0;
   /* Newton step: */ bold = b; b = b + bcov*(x'(y-mu));
   cnv = prodc(abs(b - bold) .< _bcriter) and (abs(devold - dev) < _dcriter);
endo;
/* weighted residual sum of squares: */ rss = (y-mu)'((y-mu)./mu);
/* pseudo-outcomes: z = eta + (y-mu)./mu; note that b = bcov*wx'z; */
if const gt 1 and np gt const; local k; k = seqa(const+1,1,np-const);
   /* center all covariates after the first j, then recompute bcov & b: */
   x[.,k] = x[.,k] - (mu/sumc(mu))'x[.,k];
   bcov = invpd(moment(sqrt(mu).*x,0)); b = bcov*(x'(mu.*eta + y - mu));
endif;
/* structural-robust covariance: */
   bcovr =  bcov*(x'((((y-mu))^2).*x))*bcov;
if _specrob eq 2; /* binomial covariance: */
   bcovr = bcov*(x'((y + mu.*(mu - 2*y)./n).*x))*bcov; 
endif;
if iter ge _maxit;
   "WARNING: No convergence after ";; iter;; " iterations";
endif;
if bnames[1]; format 1,0;
   if _specrob; "With inverse-information variance:"; endif;
   rreport(b,sqrt(diag(bcov)),bnames,yname,dev,rss*(-1)^(1-_fit),
           df,1+_fit*(mcc ge 3));
   if _fit; format 6,1; "The mean case count is ";; mcc;
      if mcc ge 3 and mcc lt 5;
         format 4,1; "Warning: The mean case count is less than 5, ";
         "hence the above deviance and RSS tests of fit may be invalid";
      endif;
   endif;
   if _specrob; "With structural-robust variance, assuming ";;
      if _specrob eq 2; "binomial";; else; "Poisson";; endif; " variability:";
      rreport(b,sqrt(diag(bcovr)),bnames,yname,-1,-1,df,0); endif;
   if _fit; "Total run time: ";; etstr(ethsec(t0,date)); endif;
   print; format 10,3;
endif;
retp(b,bcov,bcovr,dev,rss,df);
endp;
