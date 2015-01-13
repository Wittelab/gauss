/* This procedure does ML logistic regression with likelihood penalized
   by multiple linear constraints of the form b = z*bs, where bs is unknown.
   WARNING: THE VARIANCE ESTIMATES OF UNPENALIZED COEFFICIENTS DO NOT ACCOUNT
            FOR THE PENALTIES, AND SO ARE OVERESTIMATED.
 Inputs: x = design matrix,
         y = vector of case counts,
         n = vector of totals (if scalar, all totals will be set to this n),
         iss = index of coefficients to be modeled at second stage
              -- these need not be unique, but indices must not be repeated
                 in the same penalty, and not all df will be correct if there
                 are repeats
              -- if const = 1, the intercept will NOT be counted in the
                 indices and will NOT be modeled
              -- if is = 0, all coefficients except the intercept
                 will be modeled
         znames = global names of the 2nd-stage design matrices for the
                  linear constraints -- these are the names you have given
                  to the design matrices for each group -- the sum of the
                  number of rows of these matrices must (logically)
                  equal the length of iss
                  UNLESS _psame = 1 (see globals, below), in which case
                  the set of penalized cofficients will overlap
         NOTE: If the matrix corresponding to a name is all zeros,
               coefficients in that group will be shrunk zero, but,
               unlike logregp.g, that matrix must have as many rows as
               there are 1st-stage parameters in the group.
         t2 = second-stage residual variances -- may be
              a) a single common value to be used for all modeled groups
                 of coefficients
              b) a vector containing one value for each modeled group
              NOTE: A value of zero for a group will result in an estimate
                    of t2 being made for that group (empirical Bayes),
         rep = vector of repetion counts (0 if none),
         const = 0 if no intercept (constant),
         offset = vector of offsets (0 if no offsets),
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
         yname = scalar name of outcome (may be 0 if no printed output).
         bnames2 = vector of second-stage coefficient names
                   (set to 0 if no second-stage printed output is desired),
 Outputs: b = coefficient estimates,
          bcov = inverse-information covariance matrix,
          bs = second-stage coefficient estimates (0 if z eq 0),
          vs = variance estimates for bs (0 if z eq 0),
          t2 = t2 if input t2 gt 0, estimate of t2 if input t2 eq 0
          dev = residual deviance for model,
          devp = residual penalized deviance,
          df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _psame = 1 if the penalties all apply to the same set of coefficients,
               so iss need only list those indices once (default is zero)
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _t2init = initial value for prior variance(s) t2 (default is 1)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
      _tcriter = convergence criterion for prior variance (default is .01)
      _maxit = maximum number of iterations (default is 30)
*/
proc (8) = lrpm(x,y,n,iss,znames,t2,rep,const,offset,bnames,yname,bnames2);
declare matrix _binit = 0; declare matrix _t2init = 1;
declare matrix _bcriter = .001; declare matrix _dcriter = .01;
declare matrix _tcriter = .01; declare matrix _maxit = 30;
declare matrix _psame = 0;
local t0,nr,np,ns,eyens,eb,df2,tozero,llsat,b,bold,bcov,lpen,
      dev,devp,devpold,undsp,t2old,df,neq0,nnot0,iter,cnv,eta,p,mu,w,inf,infi,
      w2,wz,e2,sw2,e,ih,ip,s,dfp,rss,vs,bs,mcc,bz,wvce;
local nz,ie,i,z,iz,is,it,vm,sez;
t0 = date; bnames = 0$+bnames; yname = 0$+yname; bnames2 = 0$+bnames2;
if bnames[1] and rows(bnames) ne cols(x);
   "INPUT ERROR: No. 1st-stage names not equal to no. of regressors."; end;
endif;
nz = rows(znames); tozero = zeros(nz,1);
if sumc(rep); y = y.*rep; n = n.*rep; endif;
/* Delete empty covariate patterns: */ neq0 = n .eq 0;
   if sumc(neq0); x = delif(x,neq0); y = delif(y,neq0); n = delif(n,neq0);
      if rows(offset) eq rows(neq0); offset = delif(offset,neq0); endif;
   endif;
nr = rows(x);
if rows(n) eq 1; n = n*ones(nr,1); endif;
if iss[1] eq 0; iss = seqa(1,1,cols(x)); endif;
if const; /* include constant: */ x = ones(nr,1)~x; iss = iss + 1; endif;
/* No. parameters: */ np = cols(x);
if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN LRPM.G";
   retp(0,0,0,0,0,0,0,0); endif;
/* No. 2nd-stage constraints: */ ns = rows(iss);
/* Degrees of freedom without penalty: */ df = nr - np;
if rows(t2) ne 1 and rows(t2) ne nz and rows(t2) ne ns;
   "INPUT ERROR: Wrong number of input prior variances"; end; endif;
if sumc(t2 .eq 0);
   eb = (t2.*ones(nz,1)) .eq 0;
   if rows(_t2init) eq nz; t2 = eb.*_t2init + t2;
      elseif rows(_t2init) eq 1; t2 = eb.*_t2init.*ones(nz,1) + t2;
      else; "INPUT ERROR: Wrong number of initial prior variances"; end;
   endif;
   else; eb = zeros(nz,1); if rows(t2) eq 1; t2 = t2*ones(ns,1); endif;
endif;
if bnames[1] or bnames2[1]; print; format /rds 3,0;
   " Proc lrpm.g: ML logistic regression with multiple penalty functions";
   "              and ";;
   if sumc(eb); "Morris prior variance estimates";;
      if sumc(eb) lt nz; " and ";; else; "."; endif;
   endif;
   if sumc(1-eb); "fixed prior variances."; endif;
   ""$+ftocv(np,1,0)$+" first-stage regressors, of which ";;
   ""$+ftocv(rows(unique(iss,1)),1,0)$+" are modeled in the second stage";
   "using "$+ftocv(ns,1,0)$+" second-stage parameters.";
   format 5,0; sumc(n);; " total count,";; 2*nr;; " fitted cells,";;
   " and ";; format 5,2; sumc(n)/(2*nr);; " average count.";
endif;
/* Saturated loglikelihood + n'ln(n) - sumc(ln(combin(n,y))): */
   llsat = y'ln(y + (y.==0)) + (n - y)'ln(n - y + (y.==n));
/* Initialize beta, deviance, convergence & undisp indicator, penalty: */
   if rows(_binit) eq np; b = _binit;
      else; /* start with WLS estimates: */
      p = (y + sumc(y)/sumc(n))./(n + 1); w = n.*p.*(1-p);
      trap 1; bcov = invpd(moment(sqrt(w).*x,0));
      if scalerr(bcov); "SINGULARITY INITIALIZING LRPM.G";
         retp(0,0,0,0,0,0,0,0); endif; trap 0;
      b = bcov*(x'(w.*(logit(p)-offset)));
   endif;
   devp = 0; cnv = 0; t2old = t2; undsp = 0; lpen = 0; e2 = zeros(ns,1);
/* Begin iterative reweighted penalized LS estimation */ iter = 0;
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictors: */ eta = offset + x*b;
   /* expected proportions: */ p = 1./(1 + exp(-eta));
   /* expected counts: */ mu = n.*p;
   devpold = devp; dev = 2*(llsat - (y'eta + n'ln(n - mu)));
   if bnames[1] or bnames2[1];
      "Deviance, penalty, & penalized deviance at iteration ";;
      format /rds 4,0; iter;; ":"; format 12,3;"   ";; dev;; lpen;; dev+lpen;
   endif;
   /* penalize the deviance: */ devp = dev + lpen;
   /* 1st-stage weight vector: */ w = mu.*(1-p);
   /* 1st-stage residual and score: */ e = y-mu; s = x'e;
   /* information matrix: */ inf = moment(sqrt(w).*x,0);
   /* temporary storage: */ bcov = inf;
   /* inverse of unpenalized information: */ trap 1; infi = invpd(inf);
      if scalerr(infi); "SINGULARITY IN LRPM.G";
         retp(0,0,0,0,0,0,0,0); endif; trap 0;
   t2old = t2; i = 0; ie = 0; lpen = 0;
   do until i ge nz; i = i + 1;
      z = varget(""$+znames[i]); if _psame; ie = 0; endif;
      if rows(t2) eq ns; it = ie + seqa(1,1,rows(z)); else; it = 1; endif;
      iz = ie + 1; ie = ie + rows(z);
      is = iss[iz:ie]; eyens = eye(rows(z));
      /* inverse of estimated unconditional covariance of b, times z: */
         w2 = invpd(infi[is,is] + t2[it].*eyens); wz = w2*z;
      if sumall(abs(z)) eq 0; ih = eyens; df2 = rows(z);
         else; ih = eyens - z*invpd(z'wz)*wz'; df2 = rows(z) - cols(z);
      endif;
      if eb[i]; /* update Morris variance estimate: */
         /* 1-step approx to 2nd-stage residual: */
            e2[iz:ie] = infi[is,.]*s + ih*b[is]; sw2 = sumall(w2);
         t2[i] =
             (ns*(e2[iz:ie]'(w2*e2[iz:ie]))/df2 - sumall(w2*infi[is,is]))/sw2;
         t2[i] = max(t2[i],0); undsp = undsp + (t2[i] le 0);
         t2[i] = t2[i] + (t2[i] le 0)*.0001;
         if undsp gt 1; "UNDERDISPERSION in lrpm.g -";;
            "- t2[";; format 2,0; i;;"] will be fixed at 0.001.";
            t2[i] = .001; t2old[i] = t2[i]; eb[i] = 0;
         endif;
      endif;
      /* prior information: */ ip = ih'(ih./t2[it]);
      /* penalize the information: */ bcov[is,is] = inf[is,is] + ip;
      /* penalize the score: */ ip = ip*b[is]; s[is] = s[is] - ip;
      /* contribute to the LL penalty: */ lpen = lpen + b[is]'ip;
   endo;
   if (bnames[1] or bnames2[1]) and sumc(eb);
      " -- estimated 2nd-stage residual variances t2:";; selif(t2,eb)'; endif;
   /* inverse 2nd derivative of penalized LL: */
      trap 1; bcov = invpd(bcov);
      if scalerr(bcov); "PENALIZED INF SINGULAR IN LRPM.G";
         retp(0,0,0,0,0,0,0,0); endif; trap 0;
   /* Newton step: */ bold = b; b = b + bcov*s;
   cnv = prodc(abs(b - bold) .lt _bcriter)
         and (abs(devpold - devp) lt _dcriter)
         and prodc(abs(t2 - t2old) lt _tcriter);
endo;
if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
/* First-stage degrees of freedom - sumall((x*bcov).*(x.*w)) equals
   tracemat(x*bcov*(x.*w)'), but takes much less storage: */
   dfp = nr - sumall((x*bcov).*(x.*w));
/* weighted sum of squared residuals: */ rss = e'(e./w);
/* initalize 2nd-stage covariance matrix, coefficients, and expected b: */
      vs = {}; bs = {}; bz = {}; sez = {};
/* temporary storage: */ bcov = infi;
i = 0; ie = 0;
do until i ge nz; i = i + 1;
   z = varget(""$+znames[i]); if _psame; ie = 0; endif;
   if rows(t2) eq ns; it = ie + seqa(1,1,rows(z)); else; it = 1; endif;
   iz = ie + 1; ie = ie + rows(z);
   is = iss[iz:ie]; eyens = eye(rows(z));
   tozero[i] = sumall(abs(z)) eq 0; infi = bcov[is,is];
   /* inverse of estimated unconditional covariance of b, times z: */
      w2 = invpd(infi + t2[it].*eyens); wz = w2*z;
   if tozero[i]; ih = eyens; df2 = rows(z);
      bz = bz|zeros(rows(z),1); sez = sez|zeros(rows(z),1);
      else; df2 = rows(z) - cols(z);
            ih = eyens - z*invpd(z'wz)*wz';
            vm = invpd(z'wz);
            bs = bs|(vm*(wz'b[is])); vs = vs|diag(vm);
            bz = bz|(z*bs[(rows(bs)-cols(z)+1):rows(bs)]);
            sez = sez|sqrt(diag(z*vm*z'));
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
if sumc(tozero) eq nz; vs = 0; bs = 0; endif;
if bnames[1] or bnames2[1]; s = ftocv(seqa(1,1,nz),1,0);
   "Prior standard deviations";;
   if rows(t2) eq nz; " by group:"; format 12,4;
      $(((s $+ ". Est.").*eb) + ((s $+ ".Fixed").*(1-eb)))';
      else;":";;
   endif;
   sqrt(t2');
   "Number of observations with nonzero weight       :";; format 10,0; nr;
   "Number of first-stage parameters                 :";; np;
   if bnames[1];
      "With standard errors from estimated posterior covariance matrix:";
      if rows(unique(iss,1)) lt np;
         "WARNING: THERE ARE UNPENALIZED COEFFICIENTS";
         " -- THEIR STANDARD ERRORS & P-VALUES ARE OVERESTIMATED";
      endif;
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
   if bnames2[1] and sumc(tozero) lt nz; i = 0; ie = 0; ip = 0;
      do until i ge nz; i = i + 1; print; if _psame; ie = 0; endif; format 2,0;
         z = varget("" $+ znames[i]); iz = ie + 1; ie = ie + rows(z);
         if not tozero[i];
            "For group ";; i;; " -- ";; df2 = (rows(z) - cols(z))*(-1)^eb[i];
            " Estimated 2nd-stage means for the modeled coefficients,";
            "               and standard errors for these estimated means:";
            rreport(bz[iz:ie],sez[iz:ie],bnames[iz:ie],yname,-1,-1,df2,0);
            is = ip + 1; ip = ip + cols(z); format 2,0;
            "For group ";; i;; " -- Estimated 2nd-stage coefficients:";
            local y2; let y2 = "betas";
            rreport(bs[is:ip],sqrt(vs[is:ip]),bnames2[is:ip],y2,-1,-1,df2,1);
         endif;
      endo;
   endif;
   if sumc(tozero); format 2,0;
      "COEFFICIENT ESTIMATES IN 2ND-STAGE GROUP";;
      if sumc(tozero) gt 1; "S ";; else; " ";; endif;
      selif(seqa(1,1,nz),tozero)';; " WERE SHRUNK TO ZERO.";
   endif;
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(b,bcov,bs,vs,t2,dev,devp,dfp);
endp;

