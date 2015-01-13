/* Test program (see also EXPREG.TST):  output file = exprisk.out reset;
let y = 8 5 22 16; let n = 106 120 98 85; _clevel = .9;
let x[4,2] = 1 0 0 0  1 1 0 1; t2 = 1000000; let bnames = tolbutam agegt55;
let yname = death; call exprisk(x,y,n,0,1,0,bnames,yname);
call expriskp(x,y,n,0,0,t2,0,1,0,bnames,yname,0);
call exprisk(x~prodr(x),y,n,0,1,0,bnames|"tol*age",yname);
call expriskp(x~prodr(x),y,n,0,0,t2,0,1,0,bnames|"TOL*AGE",yname,0); end; */
/* This procedure does ML exponential risk regression.
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
 Outputs: b = coefficient estimates,
          bcov = inverse expected-information covariance matrix,
          bcovo = inverse observed-information covariance matrix,
          dev = residual deviance for model,
          rss = residual weighted sum of squares,
          df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
      _maxit = maximum number of iterations (default is 30)
*/
proc (6) = exprisk(x,y,n,rep,const,offset,bnames,yname);
declare matrix _binit = 0; declare matrix _bcriter = .001;
declare matrix _dcriter = .01; declare matrix _maxit = 30;
local t0,neq0,nr,np,df,llsat,b,bold,bcov,
      dev,devold,iter,cnv,eta,p,mu,w,bcovo,rss,mcc;
t0 = date; bnames = 0$+bnames; yname = 0$+yname;
if sumc(rep); y = y.*rep; n = n.*rep; endif;
/* Delete empty covariate patterns: */ neq0 = n .eq 0;
   if sumc(neq0); x = delif(x,neq0); y = delif(y,neq0); n = delif(n,neq0);
      if rows(offset) eq rows(neq0); offset = delif(offset,neq0); endif;
   endif;
nr = rows(x);
if bnames[1]; print; format /rds 3,0;
   " PROC EXPRISK.G: ML exponential risk regression for "$+yname$+".";
   if rows(bnames) ne cols(x);
      "INPUT ERROR: No. names not equal to no. of regressors."; end;
   endif;
   format 5,0; sumc(n);; " total count,";; 2*nr;; " fitted cells,";;
   " and ";; format 5,2; sumc(n)/(2*nr);; " average count.";
endif;
if const eq 1; /* add constant: */ x = ones(nr,1)~x; endif;
/* No. parameters: */ np = cols(x);
if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN EXPRISK.G";
   retp(0,0,0,0,0); endif;
/* residual degrees of freedom: */ df = nr - np;
if rows(n) eq 1; n = n*ones(nr,1); endif;
/* Saturated loglikelihood + n'ln(n) - sumc(ln(combin(n,y))): */
   llsat = y'ln(y./n + (y.==0)) + (n - y)'ln(1 - y./n + (y.==n));
/* Initialize beta, deviance, convergence indicator, counter: */
   if rows(_binit) eq np; b = _binit;
      else; /* start with WLS estimates: */
      p = (y + sumc(y)/sumc(n))./(n + 1);
      w = n.*p./(1-p);
      trap 1; bcov = invpd(moment(sqrt(w).*x,0));
      if scalerr(bcov); "SINGULARITY INITIALIZING EXPRISK.G";
         retp(0,0,0,0,0); endif; trap 0;
      b = bcov*(x'(w.*(ln(p)-offset)));
   endif;
   dev = 0; cnv = 0; iter = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictors: */ eta = offset + x*b;
   /* expected proportions: */ p = exp(eta);
   /* expected counts: */ mu = n.*p;
   devold = dev; dev = 2*(llsat - (y'eta + (n - y)'ln(1 - p + (p .ge 1))));
   if bnames[1]; "Deviance at iteration ";;
      format /rds 4,0; iter;; ":";; format 12,3; dev; endif;
   /* weight vector: */ w = mu./(1-p);
   trap 1; bcov = invpd(moment(sqrt(w).*x,0));
   if scalerr(bcov); "SINGULARITY IN EXPRISK.G";
      retp(0,0,0,0,0); endif; trap 0;
   /* Gauss-Newton step: */ bold = b; b = b + bcov*(x'((y-mu)./(1-p)));
   cnv = prodc(abs(b - bold) .lt _bcriter) and (abs(devold - dev) lt _dcriter);
endo;
if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
/* weighted sum of squared residuals: */ rss = (y-mu)'((y-mu)./(mu.*(1-p)));
if const gt 1 and np gt const; local k; k = seqa(const+1,1,np-const);
   /* center all covariates after the first j, then recompute bcov & b: */
   x[.,k] = x[.,k] - (w/sumc(w))'x[.,k];
   bcov = invpd(moment(sqrt(w).*x,0)); b = bcov*(x'(w.*(eta + y./mu - 1)));
endif;
/* invert observed information: */ bcovo = invpd(x'(((1 + (p - y./n)).*w).*x));
if bnames[1]; format 10,4;
   "Number of observations with nonzero weight:";; format 10,0; nr;
   "Number of coefficients                    :";; np;
   "With standard errors from inverse of expected information matrix:";
   rreport(b,sqrt(diag(bcov)),bnames,yname,-1,-1,df,1);
   "With standard errors from inverse of observed information matrix:";
   mcc = meanc(y)|meanc(n-y);
   rreport(b,sqrt(diag(bcovo)),bnames,yname,dev,rss,df,1+(minc(mcc) ge 3));
   "The mean case and mean noncase counts are ";; format 6,1; mcc';
   if minc(mcc) lt 5;
         "The minimum of the mean case & mean noncase counts is less than 5.";
      if minc(mcc) ge 3;
         "Hence the above deviance and rss tests of fit may be invalid.";
      endif;
   endif;
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(b,bcov,bcovo,dev,rss,df);
endp;
