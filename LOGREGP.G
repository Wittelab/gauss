/* This procedure does ML logistic regression with likelihood penalized
   by a linear constraint b = z*bs, where bs is unknown.
 Inputs: x = design matrix,
         y = vector of case counts,
         n = vector of totals (if scalar, all totals will be set to this n),
         is = index of coefficients to be modelled at second stage
              -- if const = 1, the intercept will NOT be counted in the
                 indices and will NOT be modelled
              -- if is = 0, all coefficients except the intercept
                 will be modelled
         z = second-stage design matrix for linear constraint
             (set to 0 for shrinkage of coefficients to _bprior)
         t2 = second-stage residual variances
              -- set to 0 for single estimated tau2 (empirical Bayes),
         rep = vector of repetion counts (0 if none),
         const = 0 if no intercept (constant) term in model,
         offset = vector of offsets (0 if no offsets),
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
         yname = scalar name of outcome (may be 0 if no printed output).
         bnames2 = vector of second-stage coefficient names
                  (set to 0 if no second-stage printed output is desired),
 Outputs: b = coefficient estimates,
          bcov = inverse-information covariance matrix,
          bs = second-stage coefficient estimates (0 if z eq 0),
          vs = covariance matrix for bs (0 if z eq 0),
          t2 = t2 if input t2 gt 0, estimate of t2 if input t2 eq 0
          dev = residual deviance for model,
          devp = residual penalized deviance,
          df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _bprior = prior mean for b[is] when using z = 0 (default is 0)
      _lp = vector of linear-penalty coefficients (default is 0)
      _priormu = set to one to show estimated prior means for bp (default is 0)
      _clevel = confidence level expressed as a proportion (default is .95)
      _mis = 0 if no missing values, 1 if complete-case analysis
      _mcode = missing-value scalar or column vector with element
               for each column of x (default is -99999; ignored if _mis = 0)
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _t2init = initial value for prior variance t2 (default is 1)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
      _tcriter = convergence criterion for prior variance (default is .01)
      _maxit = maximum number of iterations (default is 30)
      _bc = 1 if you want Breslow-Clayton prior variance estimate, which runs
            much more slowly and requires much more memory, but is more
            accurate in small samples (default is 0, which uses Morris
            prior variance estimate)
*/
proc (8) = logregp(x,y,n,is,z,t2,rep,const,offset,bnames,yname,bnames2);
declare matrix _priormu = 0;
declare matrix _mis = 0; declare matrix _mcode = -99999;
declare matrix _binit = 0; declare matrix _t2init = 1;
declare matrix _bcriter = .001; declare matrix _dcriter = .01;
declare matrix _tcriter = .01; declare matrix _maxit = 30;
declare matrix _bprior = 0; declare matrix _lp = 0; declare matrix _bc = 0;
local t0,neq0,nr,np,ns,eyens,mat2,t2m,t2mi,eb,df2,tozero,llsat,b,bold,bcov,
      dpen,dev,devp,devpold,undsp,t2old,df,iter,cnv,eta,p,mu,w,inf,infi,w2,wz,
      e2,sw2,e,ih,ip,s,dfp,rss,vs,bs,mcc,bz,wvce;
t0 = date; bnames = 0$+bnames; yname = 0$+yname; bnames2 = 0$+bnames2;
if bnames[1] and rows(bnames) ne cols(x);
   "INPUT ERROR: No. 1st-stage names not equal to no. of regressors."; end;
   elseif bnames2[1] and rows(bnames2) ne cols(z);
      "INPUT ERROR: No. 2nd-stage names not equal to no. of regressors."; end;
endif;
if sumc(rep); y = y.*rep; n = n.*rep; endif;
/* Delete empty covariate patterns: */ neq0 = n .eq 0;
   if _mis eq 1; /* also delete records with missing regressor values: */
      neq0 = neq0 .or sumr(x .eq _mcode'); endif;
   if sumc(neq0); x = delif(x,neq0); y = delif(y,neq0);
      if rows(n) eq rows(neq0); n = delif(n,neq0); endif;
      if rows(offset) eq rows(neq0); offset = delif(offset,neq0); endif;
   endif;
nr = rows(x);
if rows(n) eq 1; n = n*ones(nr,1); endif;
if is[1] eq 0; is = seqa(1,1,cols(x)); endif;
if const; /* add constant: */ x = ones(nr,1)~x; is = is + 1; endif;
/* No. parameters: */ np = cols(x);
if rank(x'x) lt np; "FIRST-STAGE DESIGN MATRIX RANK DEFICIENT IN LOGREGP.G";
   retp(0,0,0,0,0,0,0,0); endif;
/* Degrees of freedom without penalty: */ df = nr - np;
/* No. parameters modelled in 2nd stage: */ ns = rows(is); eyens = eye(ns);
tozero = sumall(abs(z)) eq 0;
/* matrix t2: */ mat2 = cols(t2) gt 1;
eb = t2[1,1] eq 0 and ns gt 1; if eb; t2 = _t2init; endif;
if mat2; if rank(t2) lt ns; "INPUT ERROR:";;
            "2nd-stage t2 matrix must be positive definite."; end; endif;
         t2m = t2; t2mi = invpd(t2);
   elseif sumc(t2 .le 0); "INPUT ERROR:";;
          " Pre-specified 2nd-stage t2 must be positive."; end;
   else; t2m = t2.*eyens; t2mi = eyens./t2;
endif;
if tozero; ih = eyens; df2 = ns; ip = t2mi;
   else; ip = 0; df2 = ns - cols(z);
         if rank(z'z) lt cols(z);
            "SECOND-STAGE DESIGN MATRIX RANK DEFICIENT IN LOGREGP.G";
            retp(0,0,0,0,0,0,0,0); endif;
endif;
if eb and _bc; /* for Breslow-Clayton variance: */
   local xx,dxx,u,vi,xvu,vix; xx = x[.,is]*x[.,is]'; dxx = diag(xx); endif;
if bnames[1] or bnames2[1]; print; format /rds 3,0;
   "PROC LOGREGP.G: ML logistic regression with penalty function";
   "                 ";;
   if eb;
      if  _bc; "and Breslow-Clayton prior variance estimate,";
         else; "and Morris prior variance estimate,"; endif;
      else; "and fixed prior variance,";
   endif;
   "                 with "$+ftocv(np,1,0)$+" first-stage regressors.";
   ""$+ftocv(ns,1,0)$+" 1st-stage coefficients ";;
   if tozero; "to be shrunk to ";; format 5,3; _bprior';;
      if const; " (intercept not shrunk)";; endif; ".";
      else;"to be regressed on "$+ftocv(cols(z),1,0)$+" 2nd-stage regressors.";
   endif;
   format 5,0; sumc(n);; " total count,";; 2*nr;; " fitted cells,";;
   " and ";; format 5,2; sumc(n)/(2*nr);; " average count.";
endif;
/* Saturated loglikelihood + n'ln(n) - sumc(ln(combin(n,y))): */
   llsat = y'ln(y + (y.==0)) + (n - y)'ln(n - y + (y.==n));
/* Initialize beta, deviance, convergence & underdisp indicator, counter: */
   if rows(_binit) eq np; b = _binit;
      else; /* start with WLS estimates: */
      p = (y + sumc(y)/sumc(n) + .01)./(n + 1.02);
      w = n.*p.*(1-p);
      trap 1; bcov = invpd(moment(sqrt(w).*x,0));
      if scalerr(bcov); "SINGULARITY INITIALIZING LOGREGP.G";
         retp(0,0,0,0,0,0,0,0); endif; trap 0;
      b = bcov*(x'(w.*(logit(p)-offset)));
   endif;
   devp = 0; cnv = 0; t2old = t2; undsp = 0; iter = 0;
/* Begin iterative reweighted penalized LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictors: */ eta = offset + x*b;
   /* expected proportions: */ p = 1./(1 + exp(-eta));
   if sumc(p .ge 1) or sumc(p .le 0); /* step halve: */
      b = (bold + b)/2; continue; endif;
   /* expected counts: */ mu = n.*p;
   /* deviance penalty: */
      dpen = (b[is]-_bprior)'ip*(b[is]-_bprior) + sumc(_lp.*(b[is]-_bprior));
   devpold = devp; dev = 2*(llsat - (y'eta + n'ln(n - mu)));
   /* penalize the deviance: */ devp = dev + dpen;
   if bnames[1] or bnames2[1];
      "Deviance, penalty, & penalized deviance at iteration ";;
      format /rds 4,0; iter;; ":"; format 12,3; "   ";; dev;; dpen;; devp;
   endif;
   /* 1st-stage weight vector: */ w = mu.*(1-p);
   /* residual & score: */ e = y-mu; s = x'e;
   /* 1st-stage information matrix: */ inf = moment(sqrt(w).*x,0);
   /* inverse of unpenalized information: */ trap 1; infi = invpd(inf);
      if scalerr(infi); "SINGULARITY IN LOGREGP.G";
         retp(0,0,0,0,0,0,0,0); endif; trap 0;
   /* inverse of estimated unconditional covariance of b, times z: */
       w2 = invpd(infi[is,is] + t2m); wz = w2*z;
   if not tozero; /* update 2nd-stage I-H (residual projection) matrix */
      ih = eyens - z*invpd(z'wz)*wz'; endif;
   t2old = t2;
   if eb; /* update t2: */
      if _bc; /* update Breslow-Clayton variance estimate: */
         u = (x[.,is]*(ih*b[is]) + e./w);
         vi = invpd(diagrv(xx,dxx + 1/(w*t2)));
         xvu = x[.,is]'(vi*u); vix = vi*xx;
         t2 = t2 + (ns*(xvu'xvu)/df2 - t2*tracemat(vix))/sumall(vix.*vix');
         else; /* update Morris variance estimate: */
               /* 1-step approx to 2nd-stage residual: */
                  e2 = infi[is,.]*s + ih*b[is]; sw2 = sumall(w2);
               t2 = (ns*(e2'(w2*e2))/df2 - sumall(w2*infi[is,is]))/sw2;
      endif;
      t2 = max(t2,0); undsp = undsp + (t2 le 0);
      if bnames[1] or bnames2[1];
         "-- estimated second-stage residual variance t2:";; t2; endif;
      t2 = t2 + (t2 le 0)*.001;
      if undsp gt 1; "UNDERDISPERSION in logregp.g -";;
         "- t2 will be fixed at 0.0001."; t2 = .001; t2old = t2; eb = 0; endif;
      t2m = t2*eyens; t2mi = eyens/t2;
   endif;
   /* prior information: */ ip = ih'(t2mi*ih);
   /* inverse 2nd derivative of penalized LL: */
      bcov = inf; bcov[is,is] = inf[is,is] + ip; trap 1; bcov = invpd(bcov);
      if scalerr(bcov); "PENALIZED INF SINGULAR IN LOGREGP.G";
         retp(0,0,0,0,0,0,0,0); endif; trap 0;
   /* penalized score: */ s[is] = s[is] - ip*(b[is]-_bprior) - _lp;
   /* Newton step: */ bold = b; b = b + bcov*s;
   cnv = prodc(abs(b - bold) .lt _bcriter)
         and (abs(devpold - devp) lt _dcriter)
         and (abs(t2 - t2old) lt _tcriter);
endo;
if iter ge _maxit; "WARNING: No convergence after ";;
   format 4,0; iter;; " iterations"; format 9,3;
endif;
/* First-stage degrees of freedom - sumall((x*bcov).*(x.*w)) equals
   tracemat(x*bcov*(x.*w)'), but takes much less storage: */
   dfp = nr - sumall((x*bcov).*(x.*w));
/* weighted sum of squared residuals: */ rss = e'(e./w);
/* second-stage covariance matrix, coefficients, and expected b: */
   if tozero; vs = 0; bs = 0; bz = 0;
      else; vs = invpd(z'wz); bs = vs*(wz'b[is]); bz = z*bs; endif;
local cc,zt,zs;
/* final variance computations: */
   ih = eye(np); w2 = 0*ih; w2[is,is] = t2m; w2 = invpd(infi + w2);
   if np gt ns or not tozero;
      if np gt ns; zt = zeros(np,np-ns);
         zt[delif(seqa(1,1,np),sumr(seqa(1,1,np) .eq is')),.] = eye(np-ns);
         else; zt = {};
      endif;
      if not tozero;
          zs = zeros(np,cols(z)); zs[is,.] = z; else; zs = {};
      endif;
      z = zt~zs; wz = w2*z; ih = ih - z*invpd(z'wz)*wz';
   endif;
if eb; /* add variance component to account for estimation of t2: */
   cc = (df2-2)/df2; e2 = infi*(x'(y-mu)) + ih*b;
   wvce = cc*(w2*(infi*e2));
   wvce = wvce.*sqrt((sumall(w2*infi)/sumall(w2) + t2)./(diag(infi) + t2));
   bcov = infi - cc*infi*w2*ih*infi + (2/(df2-2))*wvce*wvce';
   else; bcov = infi - infi*w2*ih*infi;
endif;
if bnames[1] or bnames2[1]; format 10,4;
   if eb; "Estimated";; else; "Specified";; endif;
   if mat2; " prior covariance matrix:";; t2;
      else; " prior standard deviation(s)            :";; sqrt(t2');
   endif;
   "Number of observations with nonzero weight       :";; format 10,0; nr;
   "Number of first-stage parameters                 :";; np;
   if not tozero;
      "Number of second-stage parameters                :";; cols(z); endif;
   if bnames[1];
      "With standard errors from estimated posterior covariance matrix:";
      rreport(b,sqrt(diag(bcov)),bnames,yname,-1,rss,0,1);
   endif;
   mcc = meanc(y)|meanc(n-y);
   "The mean case and mean noncase counts are ";; format 6,1; mcc';
   format 10,2; "Residual deviance & penalized deviance : ";; dev;; devp;
   format 7,0; "Residual & penalized degrees of freedom: ";;
   df;; format 13,2; dfp; format 10,4;
   if minc(mcc) ge 4; "                               P-values: ";;
      cdfchic(dev,df);; cdfchic(devp,dfp);
      else;
      "The minimum of the mean case and mean noncase counts is less than 4.";
      "This indicates the data are too sparse for the deviance test of fit.";
   endif;
   if not tozero; print; df2 = df2*(-1)^eb;
      if _priormu;
         "Estimated 2nd-stage means for the modelled coefficients,";
         " and standard errors for these estimated means:";
         rreport(bz,sqrt(diag(z*vs*z')),bnames[is-const],yname,-1,-1,df2,0);
      endif;
      if bnames2[1]; "Estimated 2nd-stage coefficients:";
         call rreport(bs,sqrt(diag(vs)),bnames2,0$+"betas",-1,-1,df2,1);
      endif;
   endif;
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(b,bcov,bs,vs,t2,dev,devp,dfp);
endp;

