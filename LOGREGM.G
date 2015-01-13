/* This procedure does mixed-model logistic regression via penalized
   quasi-likelihood, with penalized generalized least-squares for
   estimation of random-coefficient (prior) variances.
 Inputs: x = design matrix for random coefficients, with columns partitioned
             into g groups, i.e., x = x1~...~xm; each group has its own
             unknown scale parameter (known scale differences are in d);
             set tp[1] to zero to treat group 1 coefficients as fixed,
         y = vector of case counts,
         n = vector of totals (if scalar, all totals will be set to this n),
         g = m-vector giving the number of x-columns in groups 1,...,g
             (note that sumc(d) should equal cols(x)),
         d = scaling vector for relative sizes of variances within groups;
             should have one entry for each random coefficient (if zero,
             a common variance within each group is assumed),
         tp = m-vector of prior means for the random-coefficient variances;
              tp[1] eq 0 will make group 1 the fixed coefficients
         tv = m-vector of prior variances for the random-coefficient variances
              t2, unless _lp = 1, in which case tv should contain prior
              variances of ln(t2);
              set tv[1] to zero if group 1 contains fixed coefficients,
         c = a matrix with cols(x) columns, each row representing
             a linear combination of coefficients that will be estimated
             (set to zero if none),
         rep = vector of repetion counts (0 if none),
         const = 0 if no intercept (constant) term in model,
         offset = vector of offsets (0 if no offsets),
         bnames = vector of coefficient (x column) names
                  (set to 0 if no printed output of a is desired),
         yname = scalar name of outcome (may be 0 if no printed output).
         cnames = vector of combination names (for rows of c)
                  (set to 0 if no printed output of b is desired),
 Outputs: b = coefficient estimates for x~z,
          bcov = covariance matrix for fixed prior variances
                 WARNING: bcov is NOT corrected for variance estimation;
                 if tv has elements greater than or equal to the
                 corresponding elements in tp, diag(bcov) may severely
                 underestimate the true variances,
          bc = coefficient estimates for c (0 if c eq 0),
          bcovc = covariance matrix for bc (0 if c eq 0)
                  NOTE: A variance correction for bcovc is available
                        IF the bc are first-stage coefficients
                        (see the description of the global _vcorr)
          t2 = m-vector of variance estimates for random coefficients,
          dev = residual deviance for model,
          devp = residual penalized deviance,
          dfp = estimated residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _mis = 0 if no missing values, 1 if complete-case analysis
      _mcode = missing-value scalar or column vector with element
               for each column of x (default is -99999; ignored if _mis = 0)
      _binit = vector of initial values for a (default is weighted-least
               squares estimates of fixed coefficients and zero random
               coefficients; if user-specified and constant i requested, it
               must include a value for the constant as its first element)
      _vcorr = 1 if Morris correction to bcovc for estimation of t2 is desired
               (default is 0; the Morris correction is
               appropriate only when bc corresponds to the first-stage
               coefficients of a multistage hierarchy)
      _lp = 1 penalize ln(t2/tp), 0 penalize t2-tp (default 0)
      _bcriter = convergence criterion for a (default is .001)
      _dcriter = convergence criterion for deviance (default is .01)
      _tcriter = convergence criterion for t2 (default is .001)
      _cundsp = underdispersion criterion for t2 --
                if t2[i] < _cundsp, the coefficients in group i will be
                considered underdispersed, and t2[i] will be set to tp[i]
                -- set _cundsp > maxc(tp) to force t2 = tp
      _maxit = maximum number of iterations (default is 30)
*/
proc (8) = logregm(x,y,n,g,d,tp,tv,c,rep,const,offset,bnames,yname,cnames);
declare matrix _mis = 0; declare matrix _mcode = -99999;
declare matrix _binit = 0; declare matrix _vcorr = 0;
declare matrix _bcriter = .001; declare matrix _dcriter = .01;
declare matrix _tcriter = .001; declare matrix _cundsp = .00001;
declare matrix _lp = 0; declare matrix _maxit = 30;
local t0,np,nr,ng,fe,fi,df,llsat,b,bcov,bold,lpen,dev,devp,devpold,undsp,
      t2old,neq0,iter,cnv,eta,p,mu,w,yt,inf,dfp,rss,bc,bcovc,mcc;
t0 = date; bnames = 0$+bnames; yname = 0$+yname; cnames = 0$+cnames;
/* no. coefficients: */ np = cols(x);
if bnames[1] and rows(bnames) ne np;
   "INPUT ERROR: No. names not equal to no. of regressors."; end;
   elseif cnames[1] and rows(cnames) ne rows(c);
   "INPUT ERROR: No. comb names not equal to no. combinations."; end;
endif;
if sumc(rep); y = y.*rep; n = n.*rep; endif;
/* Delete empty covariate patterns: */ neq0 = n .eq 0;
   if _mis eq 1; /* also delete records with missing regressor values: */
      neq0 = neq0 .or sumr(x .eq _mcode'); endif;
   if sumc(neq0); x = delif(x,neq0); y = delif(y,neq0);
      if rows(n) eq rows(neq0); n = delif(n,neq0); endif;
      if rows(offset) eq rows(neq0); offset = delif(offset,neq0); endif;
   endif;
nr = rows(x); fe = tp[1] eq 0;
if rows(n) eq 1; n = n*ones(nr,1); endif;
/* indices of group 1 vars: */ fi = seqa(1,1,g[1]);
if const; /* add fixed constant: */
   if cols(c) eq np; c = zeros(rows(c),1)~c;
      elseif cols(c) ne np+1;
      "INPUT ERROR: No. comb cols not equal to no. coefficients."; end;
   endif;
   x = ones(nr,1)~x; np = np+1;
   if fe; g[1] = g[1]+1; fi = fi|g[1];
      else; g = 1|g; tp = 0|tp; fe = 1; fi = 1; endif;
endif;
/* Degrees of freedom without penalty: */ df = nr - np;
/* no. groups: */ ng = rows(g);
if fe and rank(x[.,fi]) lt g[1];
   "INPUT ERROR: FIXED-EFFECTS DESIGN MATRIX RANK DEFICIENT."; end;
   elseif (tp[1] .lt 0) or sumc(tp[2:ng] .le 0);
   "INPUT ERROR: Prior variances (t2) must be positive."; end;
   elseif rows(tp) .ne ng;
   "INPUT ERROR: Rows(tp) must equal rows(tp)."; end;
   elseif rows(tv) .ne ng;
   "INPUT ERROR: Rows(tv) must equal rows(tp)."; end;
   elseif sumc(g) .ne np;
   "INPUT ERROR: sumc(g) must equal cols(x)."; end;
endif;
local eb,ig,o,xf,bf,vf,ie,t2,itv,ud,dt,dtsel,ti,xxv,i,k,s,u,wu,wux,vwuxsq,
      cdl,vt2;
eb = _cundsp lt maxc(tp); ig = seqa(1,1,ng);
/* initial groups for variance estimation: */ undsp = (1-eb)*ones(ng,1);
   if fe; undsp[1] = 1; o = zeros(g[1],1);
      xf = x[.,fi]; itv = diagm(0|(1/tv[2:ng]));
      else; o = {}; xf = 0; bf = 0; vf = 0; itv = diagm(1/tv);
   endif;
/* initial values for estimating t2: */ t2 = tp;
dt = mexpand(eye(ng-fe),g[(1+fe):ng]); dtsel = vec(dt);
if rows(d) gt 1;
   if sumc(d .le 0); "INPUT ERROR: All elements of d must be positive."; end;
      elseif rows(d) eq np-fe*g[1];
             /* multiply dt by d inverse: */ dt = dt./d; d = o|d;
      else;
      "INPUT ERROR: Rows of d do not match no. of random coefficients."; end;
   endif;
   else; d = o|((d + (d eq 0))*ones(np-fe*g[1],1));
endif;
ti = diagm(o|selif(vec(dt.*(1/t2[(1+fe):ng]')),dtsel));
/* build variance design matrix:*/
   vwuxsq = zeros(nr^2,ng);
   xxv = vwuxsq; s = cumsumc(0|g); i = fe;
   do until i ge ng; i = i + 1; k = s[i] + seqa(1,1,g[i]);
      xxv[.,i] = vec(x[.,k]*(d[k].*x[.,k]'));
   endo;
if bnames[1] or cnames[1]; print; format /rds 3,0;
"PROC LOGREGM.G: Mixed-model logistic regression via iterative least squares,";
   if eb;
      "                with penalized-";;
      if _lp; "ratio";; else; "difference";; endif;
      " RIGLS prior variance estimates.";
      else;
      "                with fixed prior variances (semi-Bayes).";
   endif;
   "                "$+ftocv(np,1,0)$+" regressors in ";;
   ""$+ftocv(ng,1,0)$+" groups.";
   if fe; "                "$+ftocv(g[1],1,0)$+" fixed coefficient";;
      if g[1] gt 1; "s";; endif; " in group 1";;
      if const; " (including intercept)";; endif; ".";
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
      trap 1; bcov = invpd(moment(sqrt(w).*x,0) + ti);
      if scalerr(bcov); "SINGULARITY INITIALIZING LOGREGM.G";
         end; endif; trap 0;
      b = bcov*(x'(w.*(logit(p)-offset)));
   endif;
   devp = 0; cnv = 0; t2old = t2; iter = 0;
/* Begin iterative reweighted penalized LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictors: */ eta = offset + x*b;
   /* expected proportions: */ p = 1./(1 + exp(-eta));
   /* expected counts: */ mu = n.*p;
   /* penalty function: */ lpen = b'ti*b;
   devpold = devp; dev = 2*(llsat - (y'eta + n'ln(n - mu)));
   if bnames[1] or cnames[1];
      "Deviance, penalty, & penalized deviance at iteration ";;
      format /rds 4,0; iter;; ":"; format 12,3; "   ";; dev;; lpen;; dev+lpen;
   endif;
   /* penalize the deviance: */ devp = dev + lpen;
   /* 1st-stage weight vector: */ w = mu.*(1-p);
   if sumc(w .le 0); break; endif;
   /* information matrix: */ inf = moment(sqrt(w).*x,0);
   /* adjusted dependent var: */ yt = eta + (y-mu)./w;
   /* conditional covariance (given fixed coefficients & t2): */
       trap 1; bcov = invpd(inf + ti); trap 0;
       if scalerr(bcov); "SINGULAR bcov IN LOGREGM.G"; end; endif;
   bold = b; b = bcov*(x'(w.*yt));
   if fe; /* compute fixed-coefficients covariance: */
      bf = b[fi]; vf = bcov[fi,fi]; endif;
   trap 1; wu = invpd(reshape(xxv*t2,nr,nr) + diagm(1/w)); trap 0;
   if scalerr(wu); /* Warn & return if singularity occurs: */
         "SINGULAR wu IN LOGREG.M."; end; endif;
   if sumc(undsp) lt ng; /* update t2 via RIGLS: */
      ie = delif(ig,undsp);
      t2old = t2; wux = wu*x; i = 0;
      do until i ge rows(ie); i = i + 1;
         k = s[ie[i]] + seqa(1,1,g[ie[i]]);
         vwuxsq[.,ie[i]] = vec(wux[.,k]*(d[k].*wux[.,k]'));
      endo;
      u = yt - xf*bf;
      if _lp; /* penalize ln(t2/tp): */
         cdl = itv[ie,ie]*ln(t2[ie]./tp[ie]);
         vt2 = invpd(vwuxsq[.,ie]'xxv[.,ie]
                     + (itv[ie,ie] - diagm(cdl))./(t2[ie]*t2[ie]'));
         t2[ie] = t2[ie]
                  + vt2*(vwuxsq[.,ie]'(vec(u*u'-diagm(1/w)+xf*vf*xf')
                         - xxv[.,ie]*t2[ie]) - cdl./t2[ie]);
         else; /* penalize t2-tp: */
         vt2 = invpd(vwuxsq[.,ie]'xxv[.,ie] + itv[ie,ie]);
         t2[ie] = max(vt2*(vwuxsq[.,ie]'vec(u*u'-diagm(1/w)+xf*vf*xf')),0)
                  + vt2*(tp[ie]./tv[ie]);
      endif;
      if sumc(t2[ie] .le _cundsp) and (bnames[1] or cnames[1]);
         "LOGREGM.G: Underdispersion -- t2 for group ";;
         format 1,0; selif(ie,t2[ie] .le _cundsp)';;
         "now fixed at prior mean."; format 9,3;
         undsp = undsp .or (t2 .le _cundsp);
         t2 = (1-undsp).*t2 + undsp.*tp;
      endif;
      ti = diagm(o|selif(vec(dt.*(1/t2[(1+fe):ng]')),dtsel));
   endif;
   cnv = prodc(abs(b - bold) .lt _bcriter)
         and prodc(abs(c*b - c*bold) .lt _bcriter)
         and (abs(devpold - devp) lt _dcriter)
         and prodc(abs(t2 - t2old) .lt _tcriter);
endo;
if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
/* Estimated residual degrees of freedom - sumall((x*bcov).*(x.*w)) equals
   tracemat(x*bcov*(x.*w)'), but takes much less storage: */
   dfp = nr - sumall((x*bcov).*(x.*w));
/* weighted sum of squared residuals: */ rss = ((y-mu)./w)'(y-mu);
   bc = c*b; bcovc = c*bcov*c';
if _vcorr and (cols(c) gt 1) and eb;
   /* correct the covariance estimate bcovc for estimation of t2: */
   local w2,df2,infi,z,e2,wih,wv;
   w2 = c*x'wux*c'; df2 = np - fe*g[1];
   trap 1; infi = invpd(c*inf*c'); trap 0;
   if scalerr(infi); /* Warn & continue if singularity occurs: */
      "SINGULAR infi IN LOGREGM.G -- Variance correction cannot be computed.";
      else; /* add variance corrections to account for estimation of t2: */
      if fe; z = c[.,fi]; wih = pinv(z'w2*z); else; z = 0; wih = 0; endif;
      e2 = bc - z*bf; wih = w2 - w2*z*wih*z'w2;
      wv = ((df2-2)/df2)*w2*infi;
      bcovc = bcovc
              + (2/(df2-2))*wv*(e2*e2'+z*vf*z')*wv' + (2/df2)*infi*wih*infi;
   endif;
endif;
if bnames[1] or cnames[1]; ig = ig[1+fe:ng]; print;
   "Number of observations with nonzero weight: ";; format 1,0; nr;
   if eb; format 10,4;
      "Estimated prior standard deviations for the random coefficients:";
      sqrt(selif(vec(dt.*t2[(1+fe):ng]'),dtsel))'.*d[(1+g[1]):np]';
   endif;
   " t and t^2 are scale factors & scale variances for group coefficients:";
   " For group no.      ";; format 1,0; ig'; format 10,4;
   "  prior mode of t  :";; sqrt(tp[ig]');
   "  prior SD of ";; if _lp; "ln(t):";; else; "t^2  :";; endif; sqrt(tv[ig]');
   if eb;
      "  underdispersion  :";; $(0$+("yes"*undsp[ig]')$+("no"*(1-undsp[ig]')));
      "  estimated t      :";;
      else; "  pre-specified t  :";;
   endif;
   sqrt(t2[ig]');
   if fe; "  (group 1 has only fixed coefficients)"; endif;
   if bnames[1];
      "With standard errors from estimated expected covariance matrix";
      "(Warning: No variance corrections for estimation of scale factors):";
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
   if cnames[1]; print; "Estimated combinations of coefficients";;
      if _vcorr; " (with Morris variance correction)";; endif; ":";
      call rreport(bc,sqrt(diag(bcovc)),cnames,yname,-1,-1,0,1);
   endif;
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
if eb; t2 = t2.*(1 - undsp); endif;
retp(b,bcov,bc,bcovc,t2,dev,devp,dfp);
endp;
