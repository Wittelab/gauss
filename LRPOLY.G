/* This procedure does polytomous ML logistic regression
 Inputs: x = design matrix.
         y = matrix of outcome counts (do not include reference level).
         n = vector of totals,
         rep = vector of repetion count (0 if none).
         const = 0 if no intercept (constant).
         offset = scalar or column vector of offsets (0 if no offsets);
                  matrix of offsets may be used if cols(offset) = cols(y),
                  in which case a separate offset will be used for each outcome.
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired).
         ynames = outcome name or vector of outcome names, 1 for each outcome
                 category - reference (noncase) category name may be
                 omitted, but if included it should be the first name
                 (yname may be 0 if no printed output is desired).
Outputs: b = matrix coefficient estimates
             (columns correspond to index outcome levels)
         bcov = inverse-information covariance matrix for vec(b),
         dev = residual deviance for model,
         df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _collaps = 1 if data rows with identical covariates should be collapsed
                 together before fitting (default is 0)
      _clevel = confidence level expressed as a proportion (default is .95)
      _mis = 0 if no missing values, 1 if complete-case analysis
      _mcode = missing-value scalar or column vector with element
               for each column of x (default -99999; ignored if _mis = 0)
      _binit = matrix of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include values for the constants in its
               first row)
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
*/
proc (4) = lrpoly(x,y,n,rep,const,offset,bnames,ynames);
declare matrix _collaps = 0; declare matrix _clevel = .95;
declare matrix _mis = 0; declare matrix _mcode = -99999;
declare matrix _binit = 0; declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .01;
local t0,neq0,nr,np,nd,df,mcc,s,mu,p,w,b,i,llsat,bold,inf,bcov,dev,devold,
      iter,cnv,eta,den,xmu,se;
t0 = date; bnames = 0$+bnames; ynames = 0$+ynames;
if sumc((n-sumr(y)) .lt 0);
   "ERROR IN LRPOLYP.G: A record has more cases than total."; end; endif;
if sumc(rep); y = y.*rep; n = n.*rep; else; n = n.*ones(rows(x),1); endif;
nd = cols(y);
/* Delete empty covariate patterns: */ neq0 = n .eq 0;
   if _mis eq 1; /* also delete records with missing regressor values: */
      neq0 = neq0 .or sumr(x .eq _mcode'); endif;
   if sumc(neq0); x = delif(x,neq0); y = delif(y,neq0); n = delif(n,neq0);
      if rows(offset) eq rows(neq0); offset = delif(offset,neq0); endif;
   endif;
if _collaps; np = seqa(nd+2,1,cols(x));
   if bnames[1]; "Collapsing over rows with identical covariates."; endif;
   if rows(offset) eq 1; x = matcol(y~n~x,np);
      else; x = matcol(y~n~x~offset,np);
            offset = x[.,cols(x)-cols(offset)+1:cols(x)];
   endif;
   y = x[.,1:nd]; n = x[.,nd+1]; x = x[.,np];
endif;
nr = rows(x); np = cols(x);
if const; /* add constant: */ x = ones(nr,1)~x; np = 1+np; endif;
if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN LRPOLY.G";
   retp(0,0,0,0); endif;
if bnames[1]; print; format /rds 3,0;
   "PROC LRPOLY.G: polytomous logistic regression";
   if rows(bnames) ne np-const;
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
   "      with ";; nd;;" outcomes and ";; np;; " regressors";;
   if const; " (including constant)."; else; "."; endif;
   mcc = meanc(y); mcc = mcc|(meanc(n)-sumc(mcc));
   "No. of each outcome : ";; rows(x)*mcc[1:nd]';
   "No. noncases & total: ";; rows(x)*(mcc[1+nd]~sumc(mcc));
   format 6,0; "No. records used: ";; nr; "No. parameters  : ";; np*nd;
   if rows(ynames) eq 1 and nd gt 1;
      ynames = (0$+strsect(ynames,1,6)$+ftocv(seqa(1,1,nd),2,0));
      elseif rows(ynames) eq nd; ynames = ynames;
      elseif rows(ynames) eq nd+1; ynames = ynames[2:nd+1];
      else; "LRPOLY.G: Wrong number of outcome names --";
      "       supply 1 name or else 1 name for each outcome category."; end;
   endif;
endif;
if rows(_binit) eq np and cols(_binit) eq nd; b = _binit;
   if bnames[1]; "User-supplied initial values."; endif;
   else; if bnames[1]; "Intializing with weighted-least squares."; endif;
   /* temporarily, s is binary denominators: */ s = (n - sumr(y)) + y;
   mu = y + (sumc(y)./sumc(s))'; p = mu./(s + 1); w = s.*p.*(1-p);
   /* now s is initial score: */ s = x'(w.*(logit(p)-offset));
   i = 0; b = zeros(np,nd);
   do until i ge nd; i = i + 1;
      trap 1; bcov = invpd(moment(sqrt(w[.,i]).*x,0));
      if scalerr(bcov); "SINGULARITY INITIALIZING LRPOLY.G";
         retp(0,0,0,0); endif; trap 0;
      b[.,i] = bcov*s[.,i];
   endo;
endif;
 /* Degrees of freedom for grouped data: */
    df = nd*(rows(y) - np);
 /* Saturated loglikelihood - sumc(ln(combin(n,y))): */
   y = y~(n - sumr(y));
   llsat = sumall(y.*ln(y + (y.==0))) - n'ln(n);
   y = y[.,1:nd];
 /* Initialize deviance, counter, convergence criterion: */
    dev = 0; iter = 0; cnv = 0;
 do until cnv or (iter ge _maxit); iter = iter + 1;
     /* linear predictors: */ eta = offset + x*b;
     /* fitted probs: */ p = exp(eta); den = 1 + sumr(p); p = p./den;
     /* fitted case counts: */ mu = n.*p;
     devold = dev; dev = 2*(llsat - (sumall(y.*eta) - n'ln(den)));
     if bnames[1]; "Deviance at iteration ";;
        format /rds 4,0; iter;; ":";; format 12,3; dev; endif;
     if iter gt 1 and devold lt dev; dev = devold;
        /* step halve: */ b = (bold + b)/2; continue;
     endif;
     xmu = mu*~x;
     /* information matrix: */ inf = bdiagc(xmu'x) - xmu'(p*~x);
     /* covariance of b: */ trap 1; bcov = invpd(inf);
        if scalerr(bcov); 
           "INFORMATION NOT POSITIVE DEFINITE IN LRPOLY.G AT ITERATION";; iter;
           retp(0,0,0,0); 
        endif; trap 0;
     /* Newton step: */ bold = b; b = b + reshape(bcov*vec(x'(y-mu)),nd,np)';
     cnv = prodc(prodc(abs(b - bold) .lt _bcriter))
           and (abs(devold - dev) lt _dcriter);
 endo;
 if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
 endif;
 if bnames[1];
   se = reshape(sqrt(diag(bcov)),nd,np)'; i = 0;
   do until i ge nd; i = i + 1;
      call rreport(b[.,i],se[.,i],bnames,ynames[i],-1,-1,df,1);
   endo;
   "The outcome categories and their mean counts are:";
    format 8,2; $(ynames|"referent")'; mcc';
   format 8,0; "Model degrees of freedom (no. fitted parameters):";; nd*np;
   "Residual degrees of freedom (cells - parameters):";; df;
   format 10,2; "Residual deviance: ";; dev;;
   if minc(mcc) ge 3; format 7,4; ", P-value = ";; cdfchic(dev,df);
      else; "."; "The minimum mean outcome-category count is less than 3.";
      "This indicates the data are too sparse for the deviance test of fit.";
    endif;
    "Total run time: ";; etstr(ethsec(t0,date));
 endif;
retp(b,bcov,dev,df);
endp;
