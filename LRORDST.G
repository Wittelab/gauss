/* This procedure does ordinal ML logistic regression using the
   stereotype model, via constrained polytomous logistic fitting.
   With fixed scores, this is equivalent to the adjacent-category logit model.
   WARNING: With estimated scores option (s = 0), convergence can be very slow
            and ordinary standard error estimates are invalid. To address this
            problem, a naive bootstrap option is provided (set nboot > 30).
 Inputs: x = design matrix,
         y = matrix of index outcome counts (do not include reference level),
         n = vector of totals,
         s = vector of scores for index outcomes (use 0 for estimated scores)
         rep = vector of repetion counts (0 if none),
         offset = scalar or vector of offsets (0 if no offsets),
         nboot = number of bootstrap replications
                 -- operates only when using estimated scores,
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
         ynames = vector of outcome names (may be 0 if no printed output).
Outputs: b = intercept estimates (one per index category)
             and slope coefficient estimates,
         bcov = inverse-information covariance matrix for vec(b),
         dev = residual deviance for model,
         df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _maxit = maximum number of iterations (default is 200)
      _bcriter = convergence criterion for b (default is .0005)
      _dcriter = convergence criterion for deviance (default is .05)
*/
proc (4) = lrordst(x,y,n,s,rep,offset,nboot,bnames,yname);
declare matrix _maxit = 200;
declare matrix _bcriter = .0005; declare matrix _dcriter = .05;
local t0,est,nr,nc,nd,nd1,ncd,np,eyend,sumy,ys,s2,df,neq0,
      llsat,b,bold,bcov,dev,devold,iter,cnv,
      eta,den,p,mu,smu,grad,score,infi,infic,infc,se;
t0 = date; bnames = 0$+bnames; nr = rows(x); nc = cols(x);
/* Indicate if scores are to be estimated: */ est = sumc(s) eq 0;
/* No. index outcome categories: */ nd = cols(y); nd1 = nd + 1;
if sumc(rep); y = y.*rep; n = n.*rep; endif;
eyend = eye(nd);
/* Marginal outcome counts: */ sumy = sumc(y);
if bnames[1]; print; format /rdn 3,0;
   "Proc lrordst.g: Ordinal logistic regression - stereotype model";
   "     with ";; nd;;" index outcomes, ";;
   nc;;" regressors, and "; format 10,4;
   if sumc(s); "     fixed outcome scores of ";; s';
      else; "     estimated outcome scores.";
   endif;
   if rows(bnames) ne cols(x);
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
endif;
 ncd = nc + nd;
 /* no. parameters: */ np = ncd + est*nd;
 /* No. empty covariate patterns: */ neq0 = n .eq 0;
 /* Degrees of freedom: */ df = nr*nd - np + est - sumc(neq0);
 /* Saturated loglikelihood - sumc(ln(combin(n,y))): */
   y = y~(n - sumr(y));
   llsat = sumall(y.*ln(y + (y .eq 0))) - n'ln(n + neq0);
   y = y[.,1:nd];
local bp,bcovp,devp,dfp,outer,outcnv,sb,scov,sdev,sdevold,sbold,sdf;
/* get unordered lr deviance: */
{ bp,bcovp,devp,dfp } = lrpoly(x,y,n,0,1,offset,0,0);
if est; /* initialize for estimation of scores using polytomous results: */
   s = meanc(bp[2:nc+1,.]./(bp[2:nc+1,1]+(bp[2:nc+1,1] .eq 0)));
   sb = vecr(bp[1,.]'~s);
   sdev = 0; sdevold = 1;
endif;
/* initialize beta at intercept: */ b = logit(sumy/sumc(n))|zeros(nc,1);
/* Outer looping is done when the scores are estimated.
   Each outer loop involves alternating between estimating
   b holding s fixed and estimating s holding b fixed. */
outer = 0; outcnv = 0;
do until outcnv or (outer gt _maxit); outer = outer + 1;
   /* Score totals: */ ys = y*s;
   /* Squared scores: */ s2 = s.*s;
   /* Initialize deviance, counter, convergence criterion: */
   dev = 0; devold = 1; iter = 0; cnv = 0;
   do until cnv or (iter ge _maxit); iter = iter + 1;
      /* linear predictors: */ eta = b[1:nd]' + (offset + x*b[nd1:ncd])*s';
      /* fitted probs: */ p = exp(eta); den = 1 + sumr(p); p = p./den;
      /* fitted case counts: */ mu = n.*p; smu = sumc(mu);
      devold = dev; dev = 2*(llsat - (sumall(y.*eta) - n'ln(den)));
      if bnames[1]; "Deviance at inner iteration ";;
         format /rdn 4,0; iter;;
         if est; ", outer iteration ";; outer;; endif;
         ":";; format 11,3; dev;
      endif;
      /* total score: */ grad = ((sumy - smu)|(x'(ys - mu*s)));
      /* information for intercepts: */
                   infi = diagrv(eyend,smu) - mu'p;
      /* coinformation for intercepts and covariate coefficients: */
                   infic = (x'(mu.*(s' - p*s)));
      /* information for covariate coefficients: */
                   infc = x'((mu*s2 - (mu*s).*(p*s)).*x);
      bcov = invpd((infi|infic)~(infic'|infc));
      /* Newton step: */ bold = b; b = b + bcov*grad;
      cnv = prodc(prodc(abs(b - bold) .lt _bcriter))
            and (abs(devold - dev) lt _dcriter);
   endo;
   if not est; /* scores are fixed, so leave outer loop: */ break; endif;
   /* otherwise, estimate scores using polytomous lr: */
   sbold = vec(sb); sdevold = sdev;
   { sb,scov,sdev,sdf } = lrpoly(offset+x*b[nd1:ncd],y,n,0,1,0,0,0);
   s = sb[2,.]';
     outcnv = prodc(prodc(abs(vec(sb) - sbold) .lt _bcriter))
           and (abs(sdevold - sdev) lt _dcriter);
 endo;
 if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
 endif;
 if est; s = s[2:nd]; b = b|s; endif;
 if bnames[1]; /* create names for intercepts, & SEs: */
    bnames = (0$+"y = "$+ftocv(seqa(1,1,nd),2,0))|bnames;
    se = sqrt(diag(bcov));
   if est; /* add names & SEs for scores: */
      bnames = bnames|(0$+"s="$+ftocv(s,5,2));
      scov = reshape(diag(scov),2,nd);
      se = se|sqrt(scov[2:nd,2]);
   endif;
   rreport(b,se,bnames,yname,dev,-1,df,1);
   if (meanc(n) lt 5) and (df gt 0);
      "The mean denominator value (MDV) is ";; meanc(n);
      "The above deviance test of fit may not be valid if MDV < 5.";
   endif;
   "Deviance statistic against unordered model:";
   format 10,3; dev-devp;;
   format 8,0; ", df = ";; df-dfp;;
   if dev-devp gt 0;
      format 6,4;", p = ";; cdfchic(dev-devp,df-dfp);
      else; print;
   endif;
   "Total main run time: ";; etstr(ethsec(t0,date));
endif;
if est and nboot; t0 = date;
   if bnames[1];
      "Monte Carlo estimation of coefficient standard errors - ";
   endif;
   local iboot,rej,bb,bbs,yboot,i,t,bbcov,bbdev,bbdf,j,ibs;
 iboot = 0; rej = 0;
   bb = zeros(np-1,nboot);  bbs = zeros(nc*(nd-1),nboot);
   yboot = zeros(nd,nr);
   do until iboot ge nboot; iboot = iboot + 1;
      i = 0;
      do until i ge nr; i = i + 1;
         t = rndmultn((1-sumc(p[i,.]'))|p[i,.]',n[i]);
         yboot[.,i] = t[2:nd1];
      endo;
      if prodc(sumr(yboot)) eq 0; /* reject the sample: */
         iboot = iboot-1; rej = rej + 1; continue; endif;
      { bb[.,iboot],bbcov,bbdev,bbdf } = lrordst(x,yboot',n,0,0,0,0,0,0);
      bbs[.,iboot] = bb[ncd+1:np-1,iboot].*.bb[nd1:ncd,iboot];
   endo;
   bb = bb'; bbs = bbs'; bbcov = moment(bb-meanc(bb)',0)/(nboot-1);
   if bnames[1];
      se = sqrt(diag(bbcov));
      " estimates from ";; format 4,0; nboot;;" trials";;
      " (excludes ";; rej;; " rejected trials):";
      rreport(meanc(bb),se,bnames,yname,-1,-1,df,1);
      call mclimits(bb[.,nd1:np-1],bnames[nd1:np-1]);
      print; "For odds ratios:";
      call mclimits(exp(bb[.,nd1:ncd]),bnames[nd1:ncd]);
      j = 0;
      /* report results for b.*.s: */
      do until j ge nd-1; j = j + 1;
        /* indices of b*s[j] in bbs: */ ibs = seqa((j-1)*nc+1,1,nc);
        format 5,2; print; "For covariate coefficients*score ";; s[j];;":";
        rreport
        (meanc(bbs[.,ibs]),stdc(bbs[.,ibs]),bnames[nd1:ncd],yname,-1,-1,df,1);
        call mclimits(bbs[.,ibs],bnames[nd1:ncd]);
        print; "For odds ratios:";
        call mclimits(exp(bbs[.,ibs]),bnames[nd1:ncd]);
      endo;
      "Total Monte Carlo time: ";; etstr(ethsec(t0,date));
   endif;
   format 10,4; bcov = bbcov;
endif;
format 10,4;
retp(b,bcov,dev,df);
endp;
