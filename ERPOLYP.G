/* This procedure does polytomous ML exponential regression with likelihood
   penalized by a linear constraint vec(b) = z*bs, where bs is unknown.
 Inputs: x = design matrix,
         y = matrix of outcome counts,
         n = vector of person-times,
         is = indices of rows of b to be modelled at second stage
              -- if const = 1, the intercept row will NOT be counted in the
                 indices and will NOT be modelled
              -- if is = 0, all rows except the intercept row
                 will be modelled
              -- see global variable _issex (below) to model column-specific
                 elements of b
         z = second-stage design matrix for linear constraint
             (set to 0 for shrinkage of coefficients to zero)
         t2 = second-stage residual variances
              -- set to 0 for single estimated tau2 (empirical Bayes),
         rep = vector of repetion count (0 if none).
         const = 0 if no intercept (constant).
         offset = scalar or column vector of offsets (0 if no offsets);
                  matrix of offsets may be used if cols(offset) = cols(y),
                  in which case a separate offset willbe used for each outcome.
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired).
         yname = outcome name or vector of outcome names, 1 for each outcome
                 category - reference (noncase) category name may be
                 included as first name
                 (yname may be 0 if no printed output is desired).
         bnames2 = vector of second-stage coefficient names
                  (set to 0 if no second-stage printed output is desired),
Outputs: b = matrix of coefficient estimates
             (columns correspond to index outcome levels)
         bcov = inverse-information covariance matrix for vec(b),
         bs = second-stage coefficient estimates (0 if z eq 0),
         vs = covariance matrix for bs (0 if z eq 0),
         t2 = t2 if input t2 gt 0, estimate of t2 if input t2 eq 0
         dev = residual deviance for model,
         df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _binit = matrix of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include values for the constants in its
               first row)
      _t2init = initial value for prior variance t2 (default is 1)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
      _tcriter = convergence criterion for prior variance (default is .01)
      _maxit = maximum number of iterations (default is 30)
      _issex = 1 (default) if "is" (above) refers to rows of b; set to 0
               if the input "is" refers to elements of vec(b)
      WARNING: IF _issex = 0 AND const = 1, YOU MUST COUNT THE INTERCEPT
               TERMS FIRST WHEN SPECIFYING THE INDICES IN vec(b) TO MODEL
*/
proc (7) = erpolyp(x,y,n,is,z,t2,rep,const,offset,bnames,ynames,bnames2);
declare matrix _binit = 0; declare matrix _t2init = 1;
declare matrix _bcriter = .0005; declare matrix _dcriter = .01;
declare matrix _tcriter = .01; declare matrix _maxit = 30;
declare matrix _issex = 1;
local t0,neq0,nr,np,nd,sy,ns,df,eyens,eb,df2,tozero,s,mu,r,b,i,ih,ip,
      llsat,bold,bcov,bcovr,dev,lpen,devold,t2old,iter,cnv,undsp,
      eta,inf,infi,e,w2,wz,e2,dfp,j,t,vs,bs,wvce;
t0 = date; bnames = 0$+bnames; ynames = 0$+ynames; bnames2 = 0$+bnames2;
if bnames[1] and rows(bnames) ne cols(x);
   "INPUT ERROR: No. 1st-stage names not equal to no. of regressors."; end;
   elseif bnames2[1] and rows(bnames2) ne cols(z);
      "INPUT ERROR: No. 2nd-stage names not equal to no. of regressors."; end;
endif;
if sumc(rep); y = y.*rep; n = n.*rep; endif;
/* Delete empty covariate patterns: */ neq0 = n .eq 0;
   if sumc(neq0); x = delif(x,neq0); y = delif(y,neq0); n = delif(n,neq0);
      if rows(offset) eq rows(neq0); offset = delif(offset,neq0); endif;
   endif;
/* no. rows remaining: */ nr = rows(x);
/* event totals: */ sy = sumr(y);
if is[1] eq 0;
   if _issex eq 0;
      "INPUT ERROR: AT LEAST ONE OF is AND _issex MUST BE POSITIVE."; end;
   else; is = seqa(1,1,cols(x));
   endif;
endif;
if const; /* include constant: */ x = ones(rows(x),1)~x;
   if _issex; is = is + 1; endif;
endif;
/* No. parameters per index outcome: */ np = cols(x);
if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN ERPOLYP.G";
   retp(0,0,0,0,0,0,0); endif;
/* No. outcomes: */ nd = cols(y);
/* Degrees of freedom for grouped data: */ df = nd*(nr - np);
tozero = sumall(abs(z)) eq 0;
/* Actual indices and number of coefficients to be modelled: */
   if _issex;
      if maxc(is) gt np;
         "ERPOLYP.G: INDEX IN is OUT OF RANGE."; end; endif;
      if rows(z) eq rows(is) and not tozero; z = eye(nd).*.z; endif;
      is = vec(seqa(0,np,nd)' + is);
      if rows(t2) eq rows(is); t2 = ones(nd,1).*.t2;
         elseif rows(t2) ne 1; "ERPOLYP.G: Wrong size for t2."; end; endif;
      elseif maxc(is) gt np*nd;
         "ERPOLYP.G: INDEX IN is OUT OF RANGE."; end;
   endif;
   ns = rows(is);
eyens = eye(ns); eb = 0;
if sumc(t2 .le 0);
   if rows(t2) gt 1; "INPUT ERROR:";;
      " Pre-specified second-stage variances (t2) must be positive."; end;
      else; eb = 1; t2 = _t2init;
   endif;
endif;
if tozero; ih = eyens; ip = ih./t2; df2 = ns;
   else; ip = 0; df2 = ns - cols(z);
endif;
if bnames[1] or bnames2[1];
   print; format /rds 3,0;
   "PROC ERPOLYP.G: polytomous exponential regression by penalized likelihood";
   "                 ";;
   if eb; "with Morris prior variance estimate.";
           else; "with fixed prior variance."; endif;
   nd;;" outcomes and ";; np-const;; " regressors";;
   if const; " (plus constant)";; endif; ",";
   "     for a total of "$+ftocv(np*nd,1,0)$+" first-stage cofficients.";
   ""$+ftocv(ns,1,0)$+" 1st-stage coefficients ";;
   if tozero; "are to be shrunk to 0";;
      if const; " (intercepts not shrunk)";; endif; ".";
      else;"to be regressed on "$+ftocv(cols(z),1,0)$+" 2nd-stage regressor";;
           if cols(z) gt 1; "s";; endif; ".";
   endif;
   format 5,0; sumc(sy);; " total count, ";;
   (nd+1)*nr;; " fitted cells, ";;
   " and ";; format 5,2; meanc(sy)/nd;; " average count.";
   if rows(ynames) eq 1;
      ynames = (0$+strsect(ynames,1,6)$+ftocv(seqa(1,1,nd),2,0));
      elseif rows(ynames) eq nd; ynames = ynames;
      elseif rows(ynames) eq nd+1; ynames = ynames[2:nd+1];
      else; "ERPOLYP.G: Wrong number of outcome names --";
      "       supply 1 name or else 1 name for each outcome category."; end;
   endif;
   if _issex and bnames2[1] and (rows(bnames2) ne cols(z)) and not tozero;
      local tnames; tnames = zeros(rows(bnames2),nd); j = 0;
      do until j ge rows(ynames); j = j + 1; i = 0;
         do until i ge rows(bnames2); i = i + 1;
           tnames[i,j] = (0$+strsect(bnames2[i],1,4)$+strsect(ynames[j],1,4));
         endo;
      endo; bnames2 = vec(tnames);
   endif;
endif;
/* Saturated loglikelihood + sumc(sy!)): */
   llsat = sumall(y.*ln(y + (y.==0)) - y);
/* Initialize betas, deviance, counter, convergence criterion: */
   if rows(_binit) eq np and cols(_binit) eq nd; b = _binit;
      if bnames[1]; "User-supplied initial values."; endif;
      else; if bnames[1]; "Intializing with weighted-least squares."; endif;
      mu = y + n.*(sumc(y)./sumc(n))'; r = mu./n;
      /* initial score: */ s = x'(mu.*(ln(r)-offset));
      i = 0; b = zeros(np,nd);
      do until i ge nd; i = i + 1;
         trap 1; bcov = invpd(moment(sqrt(mu[.,i]).*x,0));
         if scalerr(bcov); "SINGULARITY INITIALIZING ERPOLYP.G";
            retp(0,0,0,0,0,0,0); endif; trap 0;
         b[.,i] = bcov*s[.,i];
      endo;
   endif;
   b = vec(b);
   dev = 0; iter = 0; cnv = 0; undsp = 0;
do until cnv or (iter ge _maxit); iter = iter + 1;
     /* linear predictors: */ eta = offset + x*reshape(b,nd,np)';
     /* fitted rates: */ r = exp(eta);
     /* fitted case counts: */ mu = n.*r;
     /* penalty: */ lpen = b[is]'ip*b[is];
     devold = dev; dev = 2*(llsat - sumall((y.*eta) - mu));
     if bnames[1] or bnames2[1];
        "Deviance, penalty, & penalized deviance at iteration ";;
        format /rds 4,0; iter;; ":"; format 12,3;"   ";; dev;; lpen;; dev+lpen;
     endif;
     /* penalize the deviance: */ dev = dev + lpen;
     /* 1st-stage residual: */ e = y-mu;
     /* 1st-stage information matrix: */ inf = bdiagc(x'(mu*~x));
     /* inverse of unpenalized information: */ trap 1; infi = invpd(inf);
        if scalerr(infi); "SINGULARITY IN ERPOLYP.G";
           retp(0,0,0,0,0,0,0); endif; trap 0;
     /* inverse of estimated unconditional covariance of b, times z: */
        w2 = invpd(infi[is,is] + t2.*eyens); wz = w2*z;
     if not tozero; /* update 2nd-stage I-H (residual projection) matrix */
        ih = eyens - z*invpd(z'wz)*wz'; endif;
     t2old = t2;
     if eb; /* update t2: */
        /* 1-step approx to 2nd-stage residual: */
           e2 = infi[is,.]*vec(x'e) + ih*b[is];
        t2 = max((ns*(e2'(w2*e2))/df2 - sumall(w2*infi[is,is]))/sumall(w2),0);
        undsp = undsp + (t2 le 0);
        if bnames[1];
           " Estimated 2nd-stage residual variance t2:";; t2; endif;
        t2 = t2 + (t2 le 0)*.001;
        if undsp gt 1; "UNDERDISPERSION in erpolyp.g -";;
          "- t2 will be fixed at 0.001."; t2 = .001; t2old = t2; eb = 0; endif;
      endif;
      /* prior information: */ ip = ih'(ih./t2);
      /* invert 2nd derivative of penalized LL: */
         bcov = inf; bcov[is,is] = inf[is,is] + ip;
         trap 1; bcov = invpd(bcov);
         if scalerr(bcov); "PENALIZED INF SINGULAR IN ERPOLYP.G";
            retp(0,0,0,0,0,0,0); endif; trap 0;
      /* penalized score: */ s = vec(x'e); s[is] = s[is] - ip*b[is];
      /* Newton step: */ bold = b; b = b + bcov*s;
      cnv = prodc(prodc(abs(b - bold) .lt _bcriter)) and
            (abs(devold - dev) lt _dcriter) and (abs(t2 - t2old) lt _tcriter);
endo;
/* unpenalize the deviance: */ dev = dev - lpen;
if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
/* First-stage degrees of freedom - sumall((x*bcov).*(x.*mu)) equals
   tracemat(x*bcov*(x.*mu)'), but takes much less storage: */
   dfp = 0; i = 0; s = seqa(1,1,np) - np;
   do until i ge nd; i = i + 1; s = s + np;
      dfp = dfp + sumall((x*bcov[s,s]).*(x.*mu[.,i]));
   endo;
   dfp = nd*nr - dfp;
/* second-stage covariance matrix, coefficients, and expected b: */
   if tozero; vs = 0; bs = 0;
      else; vs = invpd(z'wz); bs = vs*(wz'b[is]); endif;
bcov = infi; infi = infi[is,is];
if eb; bcov[is,is] = infi - ((df2-2)/df2)*infi*w2*ih*infi;
   /* add variance component to account for estimation of t2: */
      wvce = ((df2-2)/df2)*(w2*(infi*e2));
      wvce = wvce.*sqrt((sumall(w2*infi)/sumall(w2) + t2)./(diag(infi) + t2));
      bcov[is,is] = bcov[is,is] + (2/(df2-2))*wvce*wvce';
   else; /* fixed-tau2 covariance: */ bcov[is,is] = infi - infi*w2*ih*infi;
endif;
/* reshape b to np rows, nd columns: */ b = reshape(b,nd,np)';
if bnames[1]; local se,mcc; format 10,4;
   if eb; "Estimated";; else; "Specified";; endif;
   " prior standard deviation(s)            :";; sqrt(t2');
   format 8,0;
   "Number of covariate rows with nonzero total      :";; nr;
   "Number of first-stage parameters                 :";; nd*np;
   if not tozero;
      "Number of second-stage parameters                :";; cols(z); endif;
   if bnames[1];
      "With standard errors from estimated posterior covariance matrix:";
      se = reshape(sqrt(diag(bcov)),nd,np)'; i = 0;
      do until i ge nd; i = i + 1;
         call rreport(b[.,i],se[.,i],bnames,ynames[i],-1,-1,0,1);
      endo;
   endif;
   "The outcome categories and their mean counts are:";
    mcc = meanc(y); format 8,2; $ynames'; mcc';
   "Number of model parameters :";; format 8,0; nd*np;
   "Number of fitted cells     :";; nr*nd;
   format 10,2; "Residual deviance & penalized deviance : ";; dev;; dev+lpen;
   format 7,0; "Residual & penalized degrees of freedom: ";;
   df;; format 13,2; dfp;; format 10,4;
   if minc(mcc) ge 4; print; "                               P-values: ";;
      cdfchic(dev,df);; cdfchic(dev+lpen,dfp);
      else;
      "."; "The minimum of the mean outcome-category counts is less than 4.";
      "This indicates the data are too sparse for the deviance test of fit.";
   endif;
   if bnames2[1] and not tozero; df2 = df2*(-1)^eb;
      print; "Estimated 2nd-stage coefficients:";
      local y2; let y2 = "betas";
      call rreport(bs,sqrt(diag(vs)),bnames2,y2,-1,-1,df2,1);
   endif;
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(b,bcov,bs,vs,t2,dev,dfp);
endp;
