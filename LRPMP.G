/* This procedure does polytomous ML logistic regression with likelihood
   penalized by MULTIPLE linear constraints of the form vec(b) = z*bs,
   where bs is unknown and the constraints are independent.
   WARNING: THE VARIANCE ESTIMATES OF UNPENALIZED COEFFICIENTS DO NOT ACCOUNT
            FOR THE PENALTIES, AND SO ARE OVERESTIMATED.
 Inputs: x = design matrix,
         y = matrix of outcome counts (do not include reference level),
         n = vector of totals,
         iss = indices of rows of b to be modeled at second stage
              -- these need not be unique, but indices must not be repeated
                 in the same penalty, and not all df will be correct if there
                 are repeats
              -- if const = 1, the intercept row will NOT be counted in the
                 indices and will NOT be modeled
              -- if is = 0, all rows except the intercept row
                 will be modeled
              -- see global variable _issex (below) to model column-specific
                 elements of b
         znames = global names of the 2nd-stage design matrices for the
                  linear constraints -- these are the names you have given
                  to the design matrices for each group -- the sum of the
                  number of rows of these matrices must (logically)
                  equal the length of iss
         NOTE: If the matrix corresponding to a name is all zeros,
               coefficients in that group will be shrunk zero, but,
               unlike lrpolyp.g, that matrix must have as many rows as
               there are 1st-stage parameters in the group.
         t2 = second-stage residual variances -- may be
              a) a single common value to be used for all modeled groups
                 of coefficients
              b) a vector containing one value for each modeled group
              NOTE: A value of zero for a group will result in an estimate
                    of t2 being made for that group (empirical Bayes),
         rep = vector of repetion count (0 if none).
         const = 0 if no intercept (constant).
         offset = scalar or column vector of offsets (0 if no offsets);
                  matrix of offsets may be used if cols(offset) = cols(y),
                  in which case a separate offset willbe used for each outcome.
         bnames = vector of coefficient names
                  (set to 0 if no printed first-stage output is desired).
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
      _t2init = initial value for prior variance(s) t2 (default is 1)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
      _tcriter = convergence criterion for prior variance (default is .01)
      _maxit = maximum number of iterations (default is 30)
      _issex = 1 (default) if iss (above) refers to rows of b; set to 0
               if iss refers to elements of vec(b)
      WARNING: IF _issex = 0 AND const = 1, YOU MUST COUNT THE INTERCEPT
               TERMS WHEN SPECIFYING THE INDICES IN vec(b) TO MODEL
*/
proc (7) = lrpmp(x,y,n,iss,znames,t2,rep,const,offset,bnames,ynames,bnames2);
declare matrix _binit = 0; declare matrix _t2init = 1;
declare matrix _bcriter = .0005; declare matrix _dcriter = .01;
declare matrix _tcriter = .01; declare matrix _maxit = 30;
declare matrix _issex = 1;
local t0,neq0,nr,np,nd,ns,df,eyens,eb,df2,tozero,s,mu,p,w,b,i,ih,ip,
      llsat,bold,bcov,bcovr,dev,lpen,devold,t2old,iter,cnv,undsp,
      eta,den,xmu,inf,infi,e,w2,wz,e2,sw2,dfp,j,t,vs,bs,wvce;
local nz,ie,z,iz,is,vm;
t0 = date; bnames = 0$+bnames; ynames = 0$+ynames; bnames2 = 0$+bnames2;
if bnames[1] and rows(bnames) ne cols(x);
   "INPUT ERROR: No. 1st-stage names not equal to no. of regressors."; end;
endif;
nz = rows(znames); tozero = zeros(nz,1);
if rows(t2) ne 1 and rows(t2) ne nz;
   "INPUT ERROR: Wrong number of input prior variances"; end; endif;
if rows(_t2init) ne 1 and rows(_t2init) ne nz;
   "INPUT ERROR: Wrong number of initial prior variances"; end; endif;
if sumc(rep); y = y.*rep; n = n.*rep; endif;
/* Delete empty covariate patterns: */ neq0 = n .eq 0;
   if sumc(neq0); x = delif(x,neq0); y = delif(y,neq0); n = delif(n,neq0);
      if rows(offset) eq rows(neq0); offset = delif(offset,neq0); endif;
   endif;
/* no. rows remaining: */ nr = rows(x);
if iss[1] eq 0;
   if _issex eq 0;
      "INPUT ERROR: AT LEAST ONE OF is AND _issex MUST BE POSITIVE."; end;
   else; iss = seqa(1,1,cols(x));
   endif;
endif;
if const; /* include constant: */ x = ones(rows(x),1)~x;
   if _issex; iss = iss + 1; endif;
endif;
/* No. parameters per index outcome: */ np = cols(x);
if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN LRPMP.G";
   retp(0,0,0,0,0,0,0); endif;
/* No. outcomes: */ nd = cols(y);
/* Degrees of freedom for grouped data: */ df = nd*(nr - np);
/* Actual indices and number of coefficients to be modeled: */
   if _issex; iss = vec(seqa(0,np,nd)' + iss); endif; ns = rows(iss);
eb = (t2.*ones(nz,1)) .eq 0;
if rows(_t2init) eq nz; t2 = eb.*_t2init + t2;
   elseif rows(_t2init) eq 1; t2 = eb.*_t2init.*ones(nz,1) + t2;
endif;
if bnames[1] or bnames2[1]; print; format /rds 3,0;
   "Proc lrpmp.g: polytomous logistic regression, multiple penalty functions,";
   "              and ";;
   if sumc(eb); "Morris prior variance estimates";;
      if sumc(eb) lt nz; " and ";; else; "."; endif;
   endif;
   if sumc(1-eb); "fixed prior variances."; endif;
   nd;;" outcomes and ";; np-const;; " regressors";;
   if const; " (plus constant)";; endif; ",";
   "     for a total of ";; np*nd;; " first-stage cofficients, of which";
   ""$+ftocv(rows(unique(iss,1)),1,0)$+" are modeled in the second stage using";
   ""$+ftocv(ns,1,0)$+" second-stage parameters.";
   format 5,0; sumc(n);; " total count, ";;
   (nd+1)*nr;; " fitted cells, ";;
   " and ";; format 5,2; sumc(n)/((nd+1)*nr);; " average count.";
   if rows(ynames) eq 1;
      ynames = (0$+strsect(ynames,1,6)$+ftocv(seqa(1,1,nd),2,0));
      elseif rows(ynames) eq nd; ynames = ynames;
      elseif rows(ynames) eq nd+1; ynames = ynames[2:nd+1];
      else; "LRPMP.G: Wrong number of outcome names --";
      "       supply 1 name or else 1 name for each outcome category."; end;
   endif;
endif;
/* Saturated loglikelihood - sumc(ln(combin(n,y))): */
   y = y~(n - sumr(y));
   llsat = sumall(y.*ln(y + (y.==0))) - n'ln(n);
   y = y[.,1:nd];
/* Initialize betas, deviance, counter, convergence criterion,
   underdispersion indicator, penalty, 2nd-stage residual: */
   if rows(_binit) eq np and cols(_binit) eq nd; b = _binit;
      if bnames[1]; "User-supplied initial values."; endif;
      else; if bnames[1]; "Intializing with weighted-least squares."; endif;
      /* temporarily, s is binary denominators: */ s = (n - sumr(y)) + y;
      mu = y + (sumc(y)./sumc(s))'; p = mu./(s + 1); w = s.*p.*(1-p);
      /* now s is initial score: */ s = x'(w.*(logit(p)-offset));
      i = 0; b = zeros(np,nd);
      do until i ge nd; i = i + 1;
         trap 1; bcov = invpd(moment(sqrt(w[.,i]).*x,0));
         if scalerr(bcov); "SINGULARITY INITIALIZING LRPMP.G";
            retp(0,0,0,0,0,0,0); endif; trap 0;
         b[.,i] = bcov*s[.,i];
      endo;
   endif;
   b = vec(b);
   dev = 0; iter = 0; cnv = 0; undsp = 0; lpen = 0; e2 = zeros(rows(iss),1);
do until cnv or (iter ge _maxit); iter = iter + 1;
     /* linear predictors: */ eta = offset + x*reshape(b,nd,np)';
     /* fitted probs: */ p = exp(eta); den = 1 + sumr(p); p = p./den;
     /* fitted case counts: */ mu = n.*p;
     devold = dev; dev = 2*(llsat - (sumall(y.*eta) - n'ln(den)));
     if bnames[1] or bnames2[1];
        "Deviance, penalty, & penalized deviance at iteration ";;
        format /rds 4,0; iter;; ":"; format 12,3;"   ";; dev;; lpen;; dev+lpen;
     endif;
     /* penalize the deviance: */ dev = dev + lpen;
     /* 1st-stage residual & score: */ e = y-mu; s = vec(x'e);
     xmu = mu*~x;
     /* 1st-stage information matrix: */ inf = bdiagc(xmu'x) - xmu'(p*~x);
     /* temporary storage: */ bcov = inf;
     /* inverse of unpenalized information: */ trap 1; infi = invpd(inf);
        if scalerr(infi); "SINGULARITY IN LRPMP.G";
           retp(0,0,0,0,0,0,0); endif; trap 0;
     t2old = t2;
     i = 0; ie = 0; lpen = 0;
     do until i ge nz; i = i + 1;
        z = varget("" $+ znames[i]); iz = ie + 1; ie = ie + rows(z);
        is = iss[iz:ie]; eyens = eye(rows(z));
        /* inverse of estimated unconditional covariance of b, times z: */
           w2 = invpd(infi[is,is] + t2[i].*eyens); wz = w2*z;
        if sumall(abs(z)) eq 0; ih = eyens; df2 = rows(z);
            else; ih = eyens - z*invpd(z'wz)*wz'; df2 = rows(z) - cols(z);
        endif;
        if eb[i]; /* update Morris variance estimate: */
           /* 1-step approx to 2nd-stage residual: */
              e2[iz:ie] = infi[is,.]*s + ih*b[is];
              sw2 = sumall(w2);
           t2[i] =
              (ns*(e2[iz:ie]'(w2*e2[iz:ie]))/df2 - sumall(w2*infi[is,is]))/sw2;
           t2[i] = max(t2[i],0); undsp = undsp + (t2[i] le 0);
           t2[i] = t2[i] + (t2[i] le 0)*.0001;
           if undsp gt 1; "UNDERDISPERSION in lrpmp.g -";;
              "- t2[";; format 2,0; i;;"] will be fixed at 0.001.";
              t2[i] = .001; t2old[i] = t2[i]; eb[i] = 0;
            endif;
         endif;
         /* prior information: */ ip = ih'(ih./t2[i]);
         /* penalize the information: */ bcov[is,is] = inf[is,is] + ip;
         /* penalize the score: */ ip = ip*b[is]; s[is] = s[is] - ip;
         /* contribute to the LL penalty: */ lpen = lpen + b[is]'ip;
      endo;
      if (bnames[1] or bnames2[1]) and sumc(eb);
         " -- estimated 2nd-stage residual variances t2:";; selif(t2,eb)';
      endif;
      /* inverse 2nd derivative of penalized LL: */
      trap 1; bcov = invpd(bcov);
      if scalerr(bcov); "PENALIZED INF SINGULAR IN LRPMP.G";
         retp(0,0,0,0,0,0,0); endif; trap 0;
      /* Newton step: */ bold = b; b = b + bcov*s;
   cnv = prodc(abs(b - bold) .lt _bcriter) and (abs(devold - dev) lt _dcriter)
         and prodc(abs(t2 - t2old) lt _tcriter);
endo;
/* unpenalize the deviance: */ dev = dev - lpen;
if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
/* First-stage degrees of freedom - sumall((x*bcov).*(x.*w)) equals
   tracemat(x*bcov*(x.*w)'), but takes much less storage: */
   dfp = 0; i = 0; s = seqa(1,1,np) - np;
   do until i ge nd; i = i + 1; s = s + np;
      dfp = dfp + sumall((x*bcov[s,s]).*(xmu[.,s].*(1-p[.,i])));
      j = i; t = s;
      do until j ge nd; j = j + 1; t = t + np;
         dfp = dfp - 2*sumall((x*bcov[s,t]).*(xmu[.,s].*p[.,j]));
      endo;
   endo;
   dfp = nd*nr - dfp;
/* initalize 2nd-stage cov matrix & coefficients: */ vs = {}; bs = {};
/* temporary storage: */ bcov = infi;
i = 0; ie = 0;
do until i ge nz; i = i + 1;
   z = varget("" $+ znames[i]); iz = ie + 1; ie = ie + rows(z);
   is = iss[iz:ie]; eyens = eye(rows(z));
   tozero[i] = sumall(abs(z)) eq 0; infi = bcov[is,is];
   /* inverse of estimated unconditional covariance of b, times z: */
      w2 = invpd(infi + t2[i].*eyens); wz = w2*z;
   if tozero[i]; ih = eyens; df2 = rows(z);
      else; df2 = rows(z) - cols(z);
            ih = eyens - z*invpd(z'wz)*wz';
            vm = invpd(z'wz);
            bs = bs|(vm*(wz'b[is])); vs = vs|diag(vm);
   endif;
   if eb[i];
        bcov[is,is] = infi - ((df2-2)/df2)*infi*w2*ih*infi;
        /* add variance component to account for estimation of t2[i]: */
        wvce = ((df2-2)/df2)*(w2*(infi*e2[iz:ie]));
        wvce = wvce.*sqrt((sumall(w2*infi)/sw2 + t2[i])./(diag(infi) + t2[i]));
        bcov[is,is] = bcov[is,is] + (2/(df2-2))*wvce*wvce';
      else; /* fixed-tau2 covariance: */ bcov[is,is] = infi - infi*w2*ih*infi;
    endif;
endo;
/* reshape b to np rows, nd columns: */ b = reshape(b,nd,np)';
if sumc(tozero) eq nz; vs = 0; bs = 0; endif;
if bnames[1] or bnames2[1]; local se,mcc;
   s = ftocv(seqa(1,1,nz),1,0);
   "Prior standard deviations by group:"; format 12,4;
   $(((s $+ ". Est.").*eb) + ((s $+ ". Fixed").*(1-eb)))'; sqrt(t2');
   if bnames[1];
      if rows(unique(iss,1)) lt np;
         "WARNING: THERE ARE UNPENALIZED COEFFICIENTS";
         " -- THEIR STANDARD ERRORS & P-VALUES ARE OVERESTIMATED";
      endif;
      se = reshape(sqrt(diag(bcov)),nd,np)'; i = 0;
      do until i ge nd; i = i + 1;
         call rreport(b[.,i],se[.,i],bnames,ynames[i],-1,-1,0,1);
      endo;
   endif;
   mcc = meanc(y); mcc = (meanc(n)-sumc(mcc))|mcc;
   "The outcome categories and their mean counts are:";
    format 8,2; $("referent"|ynames)'; mcc';
   "Number of covariate rows with nonzero total:";; format 8,0; nr;
   "Number of first-stage parameters           :";; nd*np;
   "Number of fitted cells                     :";; nr*(nd+1);
   format 10,2; "Residual deviance & penalized deviance : ";; dev;; dev+lpen;
   format 7,0; "Residual & penalized degrees of freedom: ";;
   df;; format 13,2; dfp; format 10,4;
   if minc(mcc) ge 4; "                               P-values: ";;
      cdfchic(dev,df);; cdfchic(dev+lpen,dfp);
      else; "The minimum of the mean outcome-category counts is less than 4.";
      "This indicates the data are too sparse for the deviance test of fit.";
   endif;
   if sumc(tozero) lt nz; i = 0; ip = 0;
      do until i ge nz; i = i + 1; print;  format 2,0;
         if (not tozero[i]) and bnames2[1];
            z = varget("" $+ znames[i]); df2 = (rows(z) - cols(z))*(-1)^eb[i];
            is = ip + 1; ip = ip + cols(z); format 2,0;
            "For group ";; i;; " -- Estimated 2nd-stage coefficients:";
            call rreport(bs[is:ip],sqrt(vs[is:ip]),
                         bnames2[is:ip],""$+"betas",-1,-1,df2,1);
          endif;
      endo;
   endif;
   if sumc(tozero); format 2,0;
      print; "COEFFICIENT ESTIMATES IN 2ND-STAGE GROUP";;
      if sumc(tozero) gt 1; "S ";; else; " ";; endif;
      selif(seqa(1,1,nz),tozero)';; " WERE SHRUNK TO ZERO.";
   endif;
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(b,bcov,bs,vs,t2,dev,dfp);
endp;
