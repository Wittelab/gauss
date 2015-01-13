/* This program does conditional logistic regression for m:(n-m) matched
   sets with variable matching ratio, penalized by a linear constraint
   b = z*bs, where bs is unknown.
   Inputs: x = regressor matrix
           y = indicator of case status (1 = case, 0 = noncase)
           id = vector of ID numbers of matched sets to which subjects belong
           is = index of coefficients to be modelled at second stage
           z = second-stage design matrix for linear constraint
               (set to 0 for shrinkage of coefficients to zero)
           t2 = second-stage residual variances
                -- set to 0 for single estimated tau2 (empirical Bayes),
           reps = repetion count for sets BY SORTED ID NO. (set to 0 if none)
                  (reps[i] must refer to the count for the i'th ID no.)
           Note that rows(id) = rows(reps) < rows(x).
           offset = vector of offsets (set to zero if no offsets)
           mprin = 1 if you want description of matched sets, 0 otherwise
           bnames = vector of coefficient names
                    (set to 0 if no printed output is desired),
           yname = scalar name of outcome (may be 0 if no output)
           bnames2 = vector of second-stage coefficient names
                     (set to 0 if no second-stage printed output is desired).
   Outputs: b = coefficient estimates,
            bcov = inverse-information covariance matrix,
            bs = second-stage coefficient estimates (0 if z eq 0),
            vs = covariance matrix for bs (0 if z eq 0),
            t2 = t2 if input t2 gt 0, estimate of t2 if input t2 eq 0,
            dev = -2*penalized conditional loglikelihood for model.
Globals (may be set by user from calling program; otherwise defaults are used):
      _bprior = prior mean for b[is] when using z = 0 (default is 0)
      _priormu = set to one to show estimated prior means for bp (default is 0)
      _clevel = confidence level for confidence intervals (default .95)
      _mis = 0 if no missing values, 1 if complete-case analysis
      _mcode = missing-value scalar or column vector with element
               for each column of x (default -99999; ignored if _mis = 0)
      _binit = vector of initial values for b (default is unmatched
               weighted-least squares estimates)
      _t2init = initial value for prior variance t2 (default is 1)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
      _tcriter = convergence criterion for prior variance (default is .01)
      _maxit = maximum number of iterations (default is 30)
               -- if 0, there will be no updating of _binit
*/
proc (6) = condlrmp(x,y,id,is,z,t2,reps,offset,mprin,bnames,yname,bnames2);
declare matrix _priormu = 0;
declare matrix _mis = 0; declare matrix _mcode = -99999;
declare matrix _binit = 0; declare matrix _t2init = 1;
declare matrix _bcriter = .001; declare matrix _dcriter = .01;
declare matrix _tcriter = .01; declare matrix _maxit = 30;
declare matrix _bprior = 0; declare matrix _clevel = .95;
local t0,np,del,ind,uid,n,sn,m,i,nss,ids,ofs,xi,yi,xs,dev0;
local ns,eyens,mat2,t2m,t2mi,eb,df2,tozero,llsat,b,bold,bcov,lpen,
      dev,devold,undsp,t2old,df,neq0,iter,cnv,eta,mu,inf,infi,w2,wz,
      e2,sw2,e,ih,ip,s,dfp,rss,vs,bs,mcc,bz,wvce,ss;
t0 = date; bnames = 0$+bnames; yname = 0$+yname; bnames2 = 0$+bnames2;
/* No. parameters: */ np = cols(x);
if rows(offset) eq 1; offset = offset*ones(rows(x),1); endif;
/* vector of sorted indices: */ ind = sortind(id);
/* sorted covariates matrix: */ x = x[ind,.];
/* sorted outcome vector: */ y = y[ind];
/* sorted id nos.: */ id = id[ind];
/* sorted offsets: */ offset = offset[ind];
if _mis eq 1; /* delete records with missing regressor values: */
   del = sumr((x~y~id~offset) .eq _mcode') .gt 0;
   if sumc(del); x = delif(x,del); y = delif(y,del); id = delif(id,del);
      offset = delif(offset,del); del = unique(id,1);
      if rows(del) lt rows(uid);
         reps = selif(reps,sumr(uid .eq del')); uid = del; endif;
   endif;
endif;
/* For each matched set, compute the size & the no. of cases: */
   /* vector of unique id nos.: */ uid = unique(id,1);
   /* total no. in each matched set: */ n = counts(id,uid);
   /* no. of cases in each matched set: */ m = countwts(id,uid,y);
if rows(reps) eq 1; /* set reps to 1: */ reps = ones(rows(uid),1); endif;
/* Sparse-data unpenalized degrees of freedom: */ df = sumc(reps) - np;
if is[1] eq 0; is = seqa(1,1,cols(x)); endif;
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
         if rank(z) lt cols(z);
            "SECOND-STAGE DESIGN MATRIX RANK DEFICIENT IN LOGREGP.G";
            retp(0,0,0,0,0,0,0,0);
         endif;
endif;
if bnames[1]; print;
   "PROC CONDLRMP.G: penalized conditional logistic regression.";
   "                 ";;
   if eb;
      " with Morris prior variance estimate.";
       else; "with fixed prior variance.";
   endif;
   if rows(bnames) ne np;
      "INPUT ERROR: No. 1st-stage names not equal to no. of regressors."; end;
   endif;
   if bnames2[1] and rows(bnames2) ne cols(z);
      "INPUT ERROR: No. 2nd-stage names not equal to no. of regressors."; end;
   endif;
   print; format /rds 4,0;
   "No. of cases and controls   :";; format 7,0; sumc(reps.*(m~(n-m)))';
   "Number of matched sets              :";; format 10,0; sumc(reps);
   "Number of first-stage regressors    :";; np;
   if not tozero;
      "Number of second-stage regressors:";; cols(z); endif;
   if mprin; /* print description of the sets: */
      "Matched set ID no., no. cases, no. noncases, & rep. count:";
      uid~m~(n-m)~reps; print;
   endif;
   ""$+ftocv(ns,1,0)$+" 1st-stage coefficients ";;
   if tozero; "to be shrunk to ";; format 5,3; _bprior';
      else;"to be regressed on "$+ftocv(cols(z),1,0)$+" 2nd-stage regressors.";
   endif;
endif;
/* Eliminate sets with no case or no control: */ del = m.*(m-n) .eq 0;
   if sumc(del); n = delif(n,del); m = delif(m,del); uid = delif(uid,del);
      reps = delif(reps,del); del = 1-sumr(id .eq uid');
      x = delif(x,del); y = delif(y,del);
      id = delif(id,del); offset = delif(offset,del);
   endif;
if (_mis or sumc(del)) and bnames[1];
   "No. of cases and controls after eliminating sets";
   "with no complete case or no complete control record:";;
   sumc(reps.*(m~(n-m)))';
endif;
if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN CONDLRMP.G";
   retp(0,0,0,0,0,0); endif;
/* Construct a new design matrix xs, consisting of covariate sums of all
   possible subsets of size m from the n subjects in each matched set,
   minus case covariate sums, for subsets containing at least one control: */
   /* positions in x where each matched set starts: */
      sn = cumsumc(n) - n + 1;
   /* loop through the sets: */
      i = 0; nss = {}; ids = {}; xs = {}; ofs = {};
   do until i ge rows(n); i = i + 1;
      /* positions of set i members in x: */ ind = seqa(sn[i],1,n[i]);
      /* Covariate and outcome sums for all subsets of size m[i]: */
         xi = sumcomb(y[ind]~offset[ind]~x[ind,.],m[i]);
      /* Remove subsets with cases only: */ xi = selif(xi,xi[.,1] .lt m[i]);
      /* No. of subsets left: */ nss = nss|rows(xi);
      /* Matched-set id numbers for rows of xs: */
         ids = ids|(i*ones(rows(xi),1));
      /* Subtract off the case covariate totals, then paste xi to xs: */
         xs = xs|(xi[.,3:np+2] - y[ind]'x[ind,.]);
         ofs = ofs|(xi[.,2] - y[ind]'offset[ind]);
   endo;
/* Begin fitting of constant outcome to xs, with no intercept: */
   /* Note that saturated conditional loglikelhood = 0 */
   /* Initialize beta, deviance, converge & underdisp indicators, counter: */
      if rows(_binit) eq np; b = _binit;
         elseif sumc(reps) eq rows(reps);
         /* start with unmatched WLS estimates: */ local p,w,u;
         p = (y + sumc(y)/rows(y))/2; w = p.*(1-p);
         u = ones(rows(x),1);
         trap 1; bcov = invpd(moment(sqrt(w).*(u~x),0));
         if scalerr(bcov); "SINGULARITY INITIALIZING CONDLRMP.G";
            retp(0,0,0,0,0,0); endif; trap 0;
         b = bcov*((u~x)'(w.*(logit(p)-offset))); b = b[2:rows(b)];
         else; b = zeros(np,1);
      endif;
   dev = 0; devold = 1; cnv = 0; t2old = t2; undsp = 0; iter = 0;
   do until cnv or (iter ge _maxit and _maxit ne 0); iter = iter + 1;
      /* compute conditional deviance, score, & information at b: */
         devold = dev; { dev,s,inf } = mscore(xs,ofs,b,nss,reps,ids);
      /* penalty function: */ lpen = (b[is]-_bprior)'ip*(b[is]-_bprior);
      if bnames[1]; "Deviance & penalty at iteration ";;
         format /rds 4,0; iter;; ":";; format 12,3; dev;; lpen; endif;
      /* penalize the deviance: */ dev = dev + lpen;
      /* unpenalized information: */ trap 1; infi = invpd(inf);
        if scalerr(infi); "SINGULARITY IN CONDLRMP.G";
           retp(0,0,0,0,0,0); endif; trap 0;
      /* inverse of estimated unconditional covariance of b, times z: */
         w2 = invpd(infi[is,is] + t2m); wz = w2*z;
      if not tozero; /* update 2nd-stage I-H (residual projection) matrix */
         ih = eyens - z*invpd(z'wz)*wz'; endif;
      if eb; /* update Morris variance estimate: */ t2old = t2;
         /* 1-step approx to 2nd-stage residual: */
            e2 = infi[is,.]*s + ih*b[is]; sw2 = sumall(w2);
         t2 = (ns*(e2'(w2*e2))/df2 - sumall(w2*infi[is,is]))/sw2;
         t2 = max(t2,0); undsp = undsp + (t2 le 0);
         if bnames[1];
            "-- estimated second-stage residual variance t2:";; t2; endif;
          t2 = t2 + (t2 le 0)*.001;
         if undsp gt 1; "UNDERDISPERSION IN CONDLRMP.G ";;
            if bnames[1]; "-- t2 will be fixed at 0.001."; endif;
            t2 = .001; t2old = t2; eb = 0;
         endif;
         t2m = t2*eyens; t2mi = eyens/t2;
      endif;
      /* prior information: */ ip = ih'(t2mi*ih);
      /* inverse 2nd derivative of penalized conditional LL: */
         bcov = inf; bcov[is,is] = inf[is,is] + ip; bcov = invpd(bcov);
      if _maxit eq 0; break; endif;
      /* penalized score: */ s[is] = s[is] - ip*(b[is]-_bprior);
      /* Newton step: */ bold = b; b = b + bcov*s;
      cnv = prodc(abs(b - bold) .lt _bcriter) and
            (abs(devold - dev) lt _dcriter) and (abs(t2 - t2old) lt _tcriter);
   endo;
/* unpenalize the deviance: */ dev = dev - lpen;
if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
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
   cc = (df2-2)/df2; e2 = infi*s + ih*b;
   wvce = cc*(w2*(infi*e2));
   wvce = wvce.*sqrt((sumall(w2*infi)/sumall(w2) + t2)./(diag(infi) + t2));
   bcov = infi - cc*infi*w2*ih*infi + (2/(df2-2))*wvce*wvce';
   else; bcov = infi - infi*w2*ih*infi;
endif;
if bnames[1]; format 10,4;
   if eb; "Estimated";; else; "Specified";; endif;
   " prior standard deviation(s)            :";; sqrt(t2');
   "With standard errors from estimated posterior covariance matrix:";
   rreport(b,sqrt(diag(bcov)),bnames,yname,-1,-1,0,1);
   format 10,2; "Residual deviance & penalized deviance : ";; dev;; dev+lpen;
   { dev0,s,inf } = mscore(xs,ofs,zeros(np,1),nss,reps,ids);
   "Score statistic and p for all coefficients   :";;
   ss = s'invpd(inf)*s;
   format 11,3; ss;; format 11,5; cdfchic(ss,np);
   if not tozero; print; df2 = df2*(-1)^eb;
      if _priormu;
         "Estimated 2nd-stage means for the modelled coefficients,";
         " and standard errors for these estimated means:";
         rreport(bz,sqrt(diag(z*vs*z')),bnames[is],yname,-1,-1,df2,0);
      endif;
      if bnames2[1]; "Estimated 2nd-stage coefficients:";
         call rreport(bs,sqrt(diag(vs)),bnames2,0$+"betas",-1,-1,df2,1);
      endif;
   endif;
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(b,bcov,bs,vs,t2,dev);
endp;

