/* This proc does exact analysis of stratified 2-by-2 tables using
   the algorithm of Vollset et al., JASA 1991.
   Inputs: a = vector of exposed-case (a-cell) counts
           b = vector of unexposed case (b-cell) counts
           c = vector of exposed noncase (c-cell) counts
           d = vector of unexposed noncase (d-cell) counts
           names = names of disease, exposure, and stratification variables
                   (0 for no printed output)
   Outputs: rr[2] = median-unbiased point estimate of odds ratio
                  (value for which upper exact p = lower exact p)
            rr[1],rr[3] = exact confidence limits (default 95% mid-p)
            plf,puf = lower & upper Fisher null p values
            pm = 2-sided null mid-p value = 2*min(lower mid-p,upper mid-p).
WARNING: WATCH SCREEN FOR UNDERFLOW OR OVERFLOW WARNINGS --
         IF ONE OCCURS WHEN RUNNING THIS PROCEDURE, OUTPUT MAY BE JUNK;
         TRY RESETTING _oflow AND RERUNNING TO ELIMINATE WARNING
Globals (may be set by user from calling program; otherwise defaults are used):
      _clevel = confidence level (default is .95)
      _criter = convergence criterion (default is .0001)
      _maxit = maximum number of iterations (default is 1000)
      _midp = compute mid-P limits (default 1; set to 0 to get Fisher limits)
      _oflow = overflow control value (default is 708, based on Pentium)
      _rrt = "null" odds-ratio value to be tested (default is 1)
*/
proc (6) = exact2x2(a,b,c,d,names);
declare matrix _clevel = .95; declare matrix _criter = .0001;
declare matrix _maxit = 1000; declare matrix _midp = 1;
declare matrix _oflow = 708; declare matrix _rrt = 1;
local m1,st,ns,sw,oflow,ms,n1,n0,as,l,u,cl,cu,mm,k,j,i,r,lc,s,p,ip,
      plf,puf,pf,pm;
names = 0$+names;
if sumc((rows(b)|rows(c)|rows(d)) .ne rows(a));
  "INPUT ERROR IN EXACT2x2.G: Count vectors must all have same length."; end;
endif;
if (a'd eq 0) and (b'c eq 0);
  "EXACT2x2.G: No tables with nonzero diagonal or off-diagonal products.";
  retp(-1,-1,-1,1,1,1);
endif;
if names[1]; print;
   ("EXACT2x2.G: Exact 2x2 analysis of the association of "
   $+names[1]$+" and "$+names[2]$+".");
   "Crude table:";; format 10,0; sumc(a);; sumc(b);
   "            ";; sumc(c);; sumc(d); format 9,4;
   if rows(names) ge 3; format /rdn 9,4;
      "-- Analysis stratified by: ";; $names[3:rows(names)]';; "."; endif;
endif;
/* eliminate any stratum with a zero margin: */
   m1 = (a+b .eq 0) .or (a+c .eq 0) .or (b+d .eq 0) .or (c+d .eq 0);
   a = delif(a,m1); b = delif(b,m1); c = delif(c,m1); d = delif(d,m1);
/* sort strata on size: */ st = rev(sortind(a+b+c+d));
   a = a[st]; b = b[st]; c = c[st]; d = d[st];
if sumc(a) gt sumc(d); /* flip the tables: */ st = a; a = d; d = st; endif;
if sumc(b) gt sumc(c); /* flip the tables: */ st = b; b = c; c = st; endif;
sw = sumc(a) gt sumc(b); if sw; /* switch the columns: */
   st = a; a = b; b = st; st = c; c = d; d = st; _rrt = 1/_rrt; endif;
/* no. strata: */ ns = rows(a);
/* overflow control constant: */ oflow = -_oflow/ns;
/* marginal case totals: */ m1 = a + b; ms = sumc(m1);
/* exposure totals: */ n1 = a + c; n0 = b + d;
/* bounds of stratum-specific a-cells: */
   l = max(0,m1-n0); u = min(m1,n1); cl = 0|cumsumc(l); cu = 0|cumsumc(u);
/* a-cell total: */ as = sumc(a);
/* compute network coefficients: */
   /* intialize recursively defined coefficient matrix: */
      mm = zeros(ns+1,cu[ns+1]+1); mm[1,1] = 1;
   k = -1;
   do until k ge ns-1; k = k + 1; j = cl[k+2]; i = j - 1;
      do until i ge cu[k+2]; i = i + 1; r = max(l[k+1],i-cu[k+1])-1;
         do until r ge min(u[k+1],i-cl[k+1]); r = r + 1;
           mm[k+2,i+1] = (mm[k+2,i+1] + mm[k+1,i-r+1]*exp(oflow
           + lnfact(n1[k+1]) - lnfact(n1[k+1]-r) - lnfact(r)
           + lnfact(n0[k+1]) - lnfact(n0[k+1]-m1[k+1]+r) - lnfact(m1[k+1]-r)));
        endo;
      endo;
   endo;
   cl = cl[ns+1]; cu = cu[ns+1];
   lc = ln(mm[ns+1,(cl+1):(cu+1)]'); lc = lc - maxc(lc);
/* compute p-values for testing _rrt: */ s = seqa(cl,1,cu-cl+1);
/* null probs from 0 to as: */
   p = exp(lc + ln(_rrt)*s); p = p/sumc(p);
   /* Fisher p-values: */ ip = as - cl + 1;
      plf = sumc(p[1:ip]); puf = 1-plf+p[ip];
   /* 2-sided mid-p: */ pm = plf-p[ip]/2; pm = 2*min(pm,1-pm);
   if sw; _rrt = 1/_rrt; endif;
if names[1]; /* print p-values: */
   "Two-sided p-values for testing OR = ";; format 7,3; _rrt;
   "(from doubling the smaller of the lower & upper p) --";
   "                   Fisher:";; format 8,4;
   pf = 2*min(plf,puf); pf;; ";  mid-p:";; pm;
endif;
/* Find confidence limits: */
local rr,namerr,alpha,zedu,zedl,n,lrr,rad,chi,iter,cnv,
      m,targ,bm,pv,spv,pd,dpd,bold,bh,bl,sh,sl,sm,wh,wl;
rr = -(1|1|1); let namerr = "lower CL" "estimate" "upper CL";
/* alpha-level for CL: */ alpha = (1-_clevel)/2; alpha = (1-alpha)|.5|alpha;
zedu = a'd eq 0; zedl = b'c eq 0;
if zedu or zedl; /* add a small constant to a, find only one limit: */
           n = n1 + n0; a = (n.*a + n1.*m1./n)./(n + 1);
           b = m1 - a; c = n1 - a; d = n0 - b;
endif;
{ lrr,rad,chi } = ormh2x2(a,b,c,d,0);
lrr = ln(lrr); rad = cdfni(alpha[1])*rad;
i = 2*zedu;
do until i ge (3-2*zedl); i = i + 1;
   if _midp; m = 1/2; else; m = (3-i)/2; endif;
   targ = alpha[i];
   /* use bisection method on ln odds ratio: */
      /* upper & lower bracketing values: */
         bh = lrr + (i-2)*rad + 1; bl = bh - 2;
      /* pv = hypergeometric terms, sh,sl,sm are p-vals: */
         pv = exp(lc + bh*s); sh = (sumc(pv[1:ip]) - m*pv[ip])/sumc(pv);
         pv = exp(lc + bl*s); sl = (sumc(pv[1:ip]) - m*pv[ip])/sumc(pv);
      cnv = 0; iter = 0;
      do until cnv or iter ge _maxit; iter = iter + 1;
         if sl lt targ; /* push bl down & try again: */ bl = bl - 1;
            pv = exp(lc + bl*s); sl = (sumc(pv[1:ip]) - m*pv[ip])/sumc(pv);
            continue;
            elseif sh gt targ; /* push it up & try again: */ bh = bh + 1;
            pv = exp(lc + bh*s); sh = (sumc(pv[1:ip]) - m*pv[ip])/sumc(pv);
            continue;
         endif;
         bm = (bh + bl)/2;
         pv = exp(lc + bm*s); sm = (sumc(pv[1:ip]) - m*pv[ip])/sumc(pv);
         if sign(sm-targ) eq sign(sl-targ);
            /* push up bl: */ bl = bm; sl = sm;
            else; /* push down bh: */ bh = bm; sh = sm;
         endif;
         cnv = abs(bh-bl) lt _criter and abs(sm-targ) lt _criter/10;
      endo;
   if not cnv;
      if names[1]; "NO CONVERGENCE FOR ";;
         format 8,0; $upper(namerr[i]);;
         " AFTER ";; format 3,0; iter;; " ITERATIONS."; endif;
      else; rr[i] = exp(bm);
   endif;
endo;
if sw; sw = rr .eq 0; rr = (1-sw)./(rr[3 2 1]+sw); i = 4-i; endif;
if names[1]; /* print results: */
   if zedu or zedl; "Point estimate undefined. ";;
      if cnv; if _midp; "Mid-p ";; else; "Fisher-p ";; endif;
         format /rdn 2,0; 100*_clevel;; "% "$+namerr[i]$+":";;
         format /rds 10,4; rr[i]'; else; ".";
      endif;
      else; "Mid-p odds-ratio estimate and ";;
            if _midp; "mid-p ";; else; "Fisher-p ";; endif;
            format /rdn 2,0; 100*_clevel;; "% limits:";
            "          ";; format /rds 10,4;
            rr[2];;"(";;rr[1];;",";;rr[3];;")";
            if sumc(rr .lt 0); "(A -1 means the algorithm ";;
               "failed to find the estimate or limit)"; endif;
   endif;
endif;
retp(rr[2],rr[1],rr[3],plf,puf,pm);
endp;
