/* This is a GAUSS procedure for empirical-Bayes or semi-Bayes estimation
   of correlated person-time rates, or of regression coefficients from
   multiplicative relative-risk models (e.g., logistic, Cox, or
   exponential-Poisson regression). It is based on generalizations of
   Morris's 1983 (JASA) formulas (see Greenland, Stat Med 1993).
   Inputs: b = estimated log rates or coefficients (first-stage
               parameter estimates),
           v = estimated covariance matrix for b
               (may be vector of variances if the b are independent),
           z = prior design matrix (set to 0 for shrinkage to 0 vector),
           t2pass = prior variance for semi-Bayes estimation (may be vector);
                    set t2pass = 0 if prior variance is to be
                    estimated (empirical-Bayes estimation),
           bnames = vector of names associated with the components of b & bp;
                    set to 0 for no printed output,
           yname = name of first-stage regressand (outcome, disease);
                   may be 0 if there is no printed output,
           bnames2 = vector of names associated with the components of bs;
                     set to 0 if you do not want the procedure to print
                     estimates for the second-stage parameters.
    Outputs: bs = second-stage coefficient (hyperparameter) estimates,
             vs = estimated covariance matrix for bs,
             bp = estimated posterior log rates or coefficents,
             vp = estimated posterior covariance for coefficients
             t2 = t2pass for semi-Bayes, estimated prior variance for
                  empirical Bayes;
             dfhm = approximate df for hierarchical model, bounded by
                   cols(z) as tau2 goes to 0, rows(z) as tau2 goes up.
Globals (may be set by user from calling program; otherwise defaults are used):
      _bprior = prior means to shrink toward when z = 0 (default is 0)
      _priormu = set to one to show estimated prior means for bp (default is 0)
      _cc = 1 to use the curvature correction to calculate point estimates
            when prior variance is estimated -- default is 1 and recommended,
            but correction may be turned off (set to _cc zero) to give results
            in better agreement with PQL procedures
      _clevel = confidence level for confidence intervals (default .95)
      _bcriter = convergence criterion for b (default is .0005)
      _maxit = maximum number of iterations (default is 50)
*/
proc (6) = ebrr(b,v,z,t2pass,bnames,yname,bnames2);
declare matrix _priormu = 0; declare matrix _bprior = 0;
declare matrix _cc = 1;
declare matrix _maxit = 50; declare matrix _criter = .0005;
local np,tozero,eb,np2,vecv,df2,eyenp,e,t2,undisp,cnv,iter,w,ws,wv,wz,
      bs,vs,rsst,rms,t2old,bp,cc,hatw,inf,wvce,hatp,vp,it2,dfhm,dft;
bnames = 0$+bnames; bnames2 = 0$+bnames2;
/* no. 1st-stage parameters: */ np = rows(b);
/* if z = 0 shrink to zero: */ tozero = sumall(abs(z)) eq 0;
if not tozero and rows(z) ne rows(b);
   "EBRR.G: rows b and rows z are ";; rows(b);; rows(z); end; endif;
eb = sumall(t2pass) eq 0;
np2 = (1-tozero)*cols(z);
if np2 ge np-2*eb; "EBRR.G: Too many 2nd-stage parameters."; end; endif;
/* indicate vector v (indep. 1st-stage estimates): */ vecv = cols(v) eq 1;
/* 2nd-stage df: */ df2 = np - np2; eyenp = eye(np);
/* Begin estimation procedure: */
/* Subtract prior means: */ b = b - _bprior;
/* Initialize loop variables for estimation of prior variance: */
   if eb; t2 = 0; else; t2 = t2pass; endif;
if bnames[1]; print; "PROC EBRR.G:";;
   if eb; "Empirical-Bayes estimation for correlated outcomes";
          " - Morris estimate of prior variance.";
      else; "Semi-Bayes estimation for correlated outcomes.";
   endif;
   if rows(bnames) ne rows(b);
     "INPUT ERROR: No. 1st-stage names not equal to no. of coefficients."; end;
   endif;
   if bnames2[1] and rows(bnames2) ne cols(z);
     "INPUT ERROR: No. 2nd-stage names not equal to no. of regressors."; end;
   endif;
   ""$+ftocv(np,1,0)$+" coefficients shrunk toward ";; format /rds 7,3;
   if tozero; "prior means of    ";; _bprior';
      else;"estimated deviations of prior means from ";; _bprior';
           ""$+ftocv(np2,1,0)$+" estimated second-stage parameters.";
   endif;
   if eb; "Estimated";; else; "Pre-specified";; endif;
   " prior standard deviations:";; sqrt(t2');
endif;
trap 1; if vecv; inf = eyenp./v; else; inf = invpd(v); endif;
        if scalerr(inf); "V SINGULAR IN EBRR.G"; retp(0,0,0,0,0,0); endif;
trap 0;
/* Estimate second-stage coefficients (and prior variance if EB): */
   /* underdispersion & convergence indicators, and counter: */
      undisp = 0; cnv = 0; iter = 0;
   do until cnv or (iter ge _maxit); iter = iter + 1;
      /* Do prior regression: */
      /* weight matrix: */
      if vecv; w = 1/(v + t2); wv = w.*v; wz = w.*z;
         else; w = invpd(v + t2.*eyenp); wv = w*v; wz = w*z;
      endif;
      ws = sumall(w);
      /* invert 2nd-stage information: */
         if tozero; vs = 0; else; vs = invpd(z'wz); endif;
      /* 2nd-stage coefficient estimates: */ bs = vs*(wz'b);
      /* residual from 2nd-stage estimates: */ e = b - z*bs;
      /* total RSS: */
      if vecv; rsst = e'(w.*e); else; rsst = e'w*e; endif;
      /* residual mean square: */ rms = np*rsst/(df2*ws);
      /* Exit loop if prior variance is fixed in advance (semi-Bayes)
         or if there is underdispersion: */
         if not eb; break; endif;
         if undisp gt 1; " UNDERDISPERSION IN EBRR.G ";; break; endif;
       /* update t2: */ t2old = t2; t2 = max(rms - sumall(wv)/ws,0);
       /* count occurrences of underdispersion (t2 le 0): */
          undisp = undisp + (t2 le 0);
       if bnames[1]; /* give estimated prior variance: */
          format 4,0; "                   at iteration ";; iter;;":";;
          format 8,4; t2;
       endif;
       /* Test for convergence: */
          if undisp lt 1 and abs(t2-t2old) lt _criter; cnv = 1; endif;
   endo;
if iter ge _maxit;
    "WARNING: No convergence of EB after ";;
    format /rdn 4,0; iter;;" iterations;";; format 8,4;
endif;
/* weight for the prior expectations of the first-stage parameters: */
   if eb; /* use curvature correction for EB: */
       cc = max((df2 - 2)/df2,0); else; cc = 1; endif;
/* projection of b to prior mean z*bs: */ hatw = z*vs*wz';
/* hatp = projection of b to posterior mean bp;
   vp = estimated posterior covariance: */
   if vecv; vp = (eyenp - cc*wv.*(eyenp - hatw)).*v';
      else; vp = v - cc*wv'(eyenp - hatw)*v;
   endif;
if eb; /* correct 1st-stage variance for estimation of the prior variance: */
      if vecv; wvce = wv.*(e.*sqrt((sumc(wv)/ws + t2)./(v + t2)));
         else; wvce = wv*(e.*sqrt((sumall(wv)/ws + t2)./(diag(v) + t2)));
      endif;
      vp = vp + (2*cc/df2)*wvce*wvce';
endif;
   it2 = ((t2 .ne 0)./(t2 + (t2 .eq 0))).*eyenp;
   /* set cc to 1 if _cc = 0: */ cc = 1 - (1-cc)*_cc;
   hatp = invpd(inf + it2*eyenp)*(inf + (1-cc)*it2 + cc*it2*hatw);
/* EB-posterior expectations of 1st-stage parameters,
   adding back prior means: */ bp = hatp*b + _bprior;
/* Approx. df in hierarchical model: */ dfhm = tracemat(hatp);
if bnames[1]; /* print results: */ format /rds 7,3; print;
   "With standard errors from estimated posterior covariance matrix:";
   if sumc(diag(vp) .lt 0); "NEGATIVE VARIANCE ESTIMATE FOR";;
      local nv; nv = selif(seqa(1,1,np),diag(vp) .lt 0); $bnames[nv]';
   endif;
   rreport(bp,sqrt(max(diag(vp),0)),bnames,yname,-1,-1,0,1);
   df2 = df2*(-1)^(t2pass eq 0);
   if _priormu;
     "Estimated prior means for the above coefficients,";
     " and standard errors for these estimated means:";
     rreport(z*bs,sqrt(diag(z*vs*z')),bnames,yname,-1,-1,df2,0);
   endif;
   if bnames2[1];
      "Estimated 2nd-stage coefficients:";
      call rreport(bs,sqrt(diag(vs)),bnames2,"beta",-1,-1,df2,1);
   endif;
   /* total residual df */ dft = np - tracemat(hatw);
   "Total rss and df: ";; format 10,3; rsst;; dft;;
   if eb; " (should be approx. equal)";
      else; ", chi-squared p = :";; cdfchic(rsst,dft);
   endif;
endif;
retp(bs,vs,bp,vp,t2,dfhm);
endp;
