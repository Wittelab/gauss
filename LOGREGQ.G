/* This procedure does binomial-likelihood logistic regression penalized by a quadratic form
    in the coefficients, expressed as a normal prior.
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
         yname = scalar name of outcome (may be 0 if no printed output),
         is = indices of coefficients to penalize (const will be added),
         bpm = prior mean for b[is],
         bpv = prior covariance matrix or variance vector for b[is].
 Outputs: b = coefficient estimates (with intercept first if const = 1)
          bcov = inverse penalized-information covariance matrix estimate,
          bcovr = Huber robust ("sandwich") covariance matrix for b
                  -- this estimate does not depend on the binomial
                     variance assumption used to derive bcov
          devp = penalized residual deviance for model,
          rssp = penalized residual weighted sum of squares,
          dfp = penalized residual degrees of freedom.
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
*/
proc (6) = logregq(x,y,n,rep,const,offset,bnames,yname,is,bpm,bpv);
declare matrix _collaps = 0; declare matrix _clevel = .95;
declare matrix _mis = 0; declare matrix _mcode = -99999;
declare matrix _binit = 0; declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .01;
declare matrix _fit = 1; declare matrix _specrob = 0;
local t0,neq0,nr,mcc,np,ns,zcut,i,infp,llsat,b,bold,bcov,devp,dev,devpold,dpen,
      df,iter,cnv,eta,p,mu,w,dfp,rssp,bcovr;
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
nr = rows(x); np = cols(x); ns = rows(bpm);
if is eq 0; is = seqa(1,1,ns); endif;
if bnames[1]; print; 
  "PROC LOGREGPQ.G: PL logistic regression for "$+yname$+" with quadratic-form penalty.";
   if rows(bnames) ne np;
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
   mcc = meanc(y)|meanc(n-y);
   print; "Nos. of cases and noncases :";; format /rds 6,0; rows(x)*mcc';
   "No. records used: ";; nr; "No. parameters  : ";; (const eq 1)+np;
   "Prior      b medians, variances, odds ratios, "$+ftocv(100*_clevel,1,0)$+"% intervals";;  
   if cols(bpv) gt 1;  
      " (correlated priors):"; sep = sqrt(diag(bpv));
      else; sep = sqrt(bpv); ":";
   endif;
   zcut = invnorm(.5+_clevel/2)*(0~-1~1); i=0;
   do until i ge ns; i= i+1;
        format /lds 8,0; $bnames[is[i]];; format /rds 9,3; bpm[i];; sep[i]^2;; 
        exp(bpm[i]+zcut*sep[i]);
   endo; print;
endif;
if const eq 1; /* add constant: */ x = ones(nr,1)~x; np = 1+np; is = 1+is; endif;
b = zeros(np,1); b[is] = bpm; bpm = b;
infp = zeros(np,np);  if cols(bpv) eq 1; bpv = diagm(bpv); endif; infp[is,is] = invpd(bpv);
if rank(x'x+infp) lt np; "DESIGN+INF MATRIX RANK DEFICIENT IN LOGREGPQ.G";
   retp(0,0,0,0,0,0); endif;
/* Degrees of freedom before penalty: */ df = nr - np;
/* Saturated loglikelihood + n'ln(n) - sumc(ln(combin(n,y))): */
   llsat = y'ln(y + (y.==0)) + (n - y)'ln(n - y + (y.==n));
/* Initialize beta, deviance, convergence indicator, counter: */
   if rows(_binit) eq np; b = _binit;
      else; /* start with WLS estimates: */
      p = (y + sumc(y)/sumc(n) + .01)./(n + 1.02);
      w = n.*p.*(1-p);
      trap 1; bcov = invpd(moment(sqrt(w).*x,0)+infp);
      if scalerr(bcov); "SINGULARITY INITIALIZING LOGREGPQ.G";
         retp(0,0,0,0,0,0); endif; trap 0;
      b = bcov*(x'(w.*(logit(p)-offset)));
   endif;
devp = 0; dev = 0; cnv = 0; iter = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictors: */ eta = offset + x*b;
   /* expected proportions: */ p = 1./(1 + exp(-eta));
   if sumc(p .ge 1) or sumc(p .le 0); /* step halve: */
      b = (bold + b)/2; continue; endif;
   /* expected counts: */ mu = n.*p;
   /* deviance penalty: */  dpen = (b-bpm)'infp*(b-bpm);
   devpold = devp; dev = 2*(llsat - (y'eta + n'ln(n - mu)));
   /* penalize the deviance: */ devp = dev + dpen;
   if bnames[1];
      "Deviance, penalty, & penalized deviance at iteration "$+ftocv(iter,1,0)$+": ";;
      format 12,3; "   ";; dev;; dpen;; devp;
   endif;
   if iter gt 1 and devpold lt devp; devp = devpold;
      /* step halve: */ b = (bold + b)/2; continue;
   endif;
   /* weight vector: */ w = mu.*(1-p);
   /* covariance of b: */ trap 1; bcov = invpd(moment(sqrt(w).*x,0)+infp);
      if scalerr(bcov); "INFORMATION NOT POSITIVE DEFINITE IN LOGREGPQ.G";
         retp(0,0,0,0,0,0); endif; trap 0;
   /* Newton step: */ bold = b; b = b + bcov*((x'(y-mu)) - infp*(b-bpm));
   cnv = prodc(abs(b - bold) .lt _bcriter) and (abs(devpold - devp) lt _dcriter);
endo;
if iter ge _maxit; "WARNING: No convergence after ";;
   format 4,0; iter;; " iterations"; format 9,3;
endif;
/* weighted sum of squared residuals: */ rssp = (y-mu)'((y-mu)./w) + dpen;
/* pseudo-residual: e = (y - mu)./w; */
/* pseudo-observed logit: z = eta + e; */
/* note that b = bcov*wx'z; */
if const gt 1 and np gt const; local k; k = seqa(const+1,1,np-const);
   /* center all covariates after the first j, then recompute bcov & b: */
   x[.,k] = x[.,k] - (w/sumc(w))'x[.,k];
   bcov = invpd(moment(sqrt(w).*x,0)+bpv); b = bcov*(x'(w.*eta + y - mu));
endif;
/* First-stage degrees of freedom - sumall((x*bcov).*(x.*w)) equals
   tracemat(x*bcov*(x.*w)'), but takes much less storage: */
   dfp = nr - sumall((x*bcov).*(x.*w));
/* the squared residuals are y.*(1-p)^2 + (n-y).*p^2 = y + p.*(mu - 2*y) */
/* hence the robust ("sandwich") covariance is: */
   bcovr = bcov*(x'((y.*(1-2*p)+ mu.*p).*x))*bcov; 
if bnames[1];
   if _specrob; "With inverse-penalized information variances:"; endif;
   rreport(b,sqrt(diag(bcov)),bnames,yname,dev,rssp*(-1)^(1-_fit),dfp,1+_fit*(minc(mcc) ge 3));
   if _specrob; print; "With specification-robust variances:";
      rreport(b,sqrt(diag(bcovr)),bnames,yname,-1,-1,df,1); endif;
   if _fit; "Total run time: ";; etstr(ethsec(t0,date)); endif;
   print; format 10,3;
endif;
retp(b,bcov,bcovr,devp,rssp,dfp);
endp;

