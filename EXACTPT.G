/* This proc does exact analysis of stratified person-time rate data using
   the algorithm of Guess & Thomas, Epidemiology 1990.
   Inputs: a = vector of exposed-case (a-cell) counts
           b = vector of unexposed case (b-cell) counts
           t1 = vector of exposed person-time
           t0 = vector of unexposed person-time
           names = names of disease, exposure, and stratification variables
                   (0 for no printed output)
   Outputs: rrm = median-unbiased point estimate of rate ratio
                  (value for which upper exact p = lower exact p)
            rrl,rru = exact confidence limits (default 95% mid-p)
            plf,puf = lower & upper Fisher null p values
            pm = 2-sided null mid-p value = 2*min(lower mid-p,upper mid-p).
WARNING: WATCH SCREEN FOR UNDERFLOW OR OVERFLOW WARNINGS --
         IF ONE OCCURS WHEN RUNNING THIS PROCEDURE, OUTPUT MAY BE JUNK
Globals (may be set by user from calling program; otherwise defaults are used):
      _clevel = confidence level (default is .95)
      _criter = convergence criterion (default is .0001)
      _maxit = maximum number of iterations (default is 1000)
      _midp = compute mid-P limits (default 1; set to 0 to get Fisher limits)
      _rrt = "null" rate-ratio value to be tested (default is 1)
*/
proc (6) = exactpt(a,b,t1,t0,names);
declare matrix _clevel = .95; declare matrix _criter = .0001;
declare matrix _maxit = 1000; declare matrix _midp = 1;
declare matrix _rrt = 1;
local ns,m1,ms,sw,as,ptr,v,mm,k,i,r,lc,s,p,plf,puf,pf,pm,rrm,rrl,rru;
names = 0$+names;
m1 = (a+b).*t1.*t0 .gt 0;
if sumc(m1) eq 0;
   "EXACTPT.G: No tables with nonzero diagonal or off-diagonal products.";
  retp(-1,-1,-1,1,1,1);
endif;
/* eliminate strata with no cases or zero denominators: */
   a = selif(a,m1); b = selif(b,m1); t1 = selif(t1,m1); t0 = selif(t0,m1);
/* no. strata: */ ns = rows(a);
if names[1]; print;
   ("EXACTPT.G: Exact analysis of the association of "
   $+names[1]$+" and "$+names[2]$+".");
   "Crude table:";; format 10,0; sumc(a);; sumc(b);
   "            ";; format 10,2; sumc(t1);; sumc(t0); format 9,4;
if sumc(a) and sumc(b); 
    "Crude rate ratio: ";; sumc(a)*sumc(t0)/(sumc(b)*sumc(t1)); 
endif;
   if rows(names) ge 3; format /rdn 9,4;
      "-- Analysis stratified by: ";; $names[3:rows(names)]';; "."; endif;
endif;
/* marginal case totals: */ m1 = a + b; ms = sumc(m1);
/* indicate if a is closer to the case total than to 0: */
   sw = 2*sumc(a) gt ms;
if sw; /* reverse columns: */ _rrt = 1/_rrt;
   ptr = a; a = b; b = ptr; ptr = t1; t1 = t0; t0 = ptr;
endif;
/* a-cell total: */ as = sumc(a);
/* person-time ratio: */ ptr = t1./t0;
/* compute network coefficients: */
   /* for bounds: */ v = min(as,cumsumc(m1)); v = 0|v;
   /* intialize recursively defined coefficient matrix: */
      mm = zeros(ns+1,as+1); mm[1,1] = 1;
   k = -1;
   do until k ge ns-1; k = k + 1; i = -1;
      do until i ge as; i = i + 1; r = max(0,i-v[k+1])-1;
         do until r ge min(i,m1[k+1]); r = r + 1;
            mm[k+2,i+1] =
               (mm[k+2,i+1] + mm[k+1,i-r+1]*combin(m1[k+1],r)*(ptr[k+1]^r));
         endo;
       endo;
    endo;
    lc = ln(mm[ns+1,.]');
/* compute p-values for testing _rrt: */ s = seqa(0,1,as+1);
/* null probs from 0 to as: */
   p = exp(-m1'ln(_rrt*ptr+1) + lc + s*ln(_rrt));
   /* Fisher p-values: */
      if sw; puf = sumc(p); plf = 1-puf+p[as+1];
           else; plf = sumc(p); puf = 1-plf+p[as+1]; endif;
   /* 2-sided p: */ pf = 2*min(plf,puf);
      pm = plf-p[as+1]/2; pm = 2*min(pm,1-pm);
if names[1]; /* print p-values: */
   "Two-sided p-values for testing RR = ";; format 5,2; _rrt;
   "(from doubling the smaller of the lower & upper p) --";
   "                   Fisher:";; format 8,4; pf;;
   ";  mid-p:";; pm;
endif;
/* Find confidence limits: */
local rr,namerr,alpha,zed,rdmh,dse,lrr,lse,chi,iter,cnv,
      m,st,bm,pv,pd,dpd,bold,bh,bl,sh,sl,sm,wh,wl;
rr = -(1|1|1); let namerr = "lower CL" "estimate" "upper CL";
/* alpha-level for CL: */ alpha = (1-_clevel)/2;
zed = sumc(a) eq 0;
if zed; /* add a small constant to a for MH, find only upper limit: */
   a = t1.*m1./((t1+t0).*(m1+1)); b = m1 - a; i = 2;
   else; i = 0;
endif;
{ rdmh,dse,lrr,lse,chi } = mhrate(a,b,t1,t0,0);
do until i ge 3; i = i + 1;
   if _midp; m = 1/2; else; m = (3-i)/2; endif;
   st = (i eq 3)*alpha + (i eq 2)/2 + (i eq 1)*(1-alpha);
   /* use bisection method on ln rate ratio: */
      /* upper & lower bracketing values: */
         bh = lrr + (i-2)*cdfni(1-alpha)*lse + 1; bl = bh - 2;
      /* pv = vector of a-cell probs 0 to as; sh,sl,sm are p-vals: */
         pv = exp(-m1'ln(exp(bh)*ptr+1) + lc + s*bh);
         sh = sumc(pv) - m*pv[as+1];
         pv = exp(-m1'ln(exp(bl)*ptr+1) + lc + s*bl);
         sl = sumc(pv) - m*pv[as+1];
      cnv = 0; iter = 0;
      do until cnv or iter ge _maxit; iter = iter + 1;
         if sl lt st; /* push it down and try again: */ bl = bl - 1;
               pv = exp(-m1'ln(exp(bl)*ptr+1) + lc + s*bl);
               sl = sumc(pv) - m*pv[as+1]; continue;
            elseif sh gt st; /* push it up & try again: */ bh = bh + 1;
               pv = exp(-m1'ln(exp(bh)*ptr+1) + lc + s*bh);
               sh = sumc(pv) - m*pv[as+1]; continue;
         endif;
         bm = (bh + bl)/2;
         pv = exp(-m1'ln(exp(bm)*ptr+1) + lc + s*bm);
         sm = sumc(pv) - m*pv[as+1];
         if sign(sm-st) eq sign(sl-st);
            /* push up bl: */ bl = bm; sl = sm;
            else; /* push down bh: */ bh = bm; sh = sm;
         endif;
         cnv = abs(bh-bl) lt _criter and abs(sm-st) lt _criter/10;
      endo;
   if not cnv;
      if names[1]; "NO CONVERGENCE FOR ";;
         format 8,0; $upper(namerr[i+(4-2*i)*sw]);;
         " AFTER ";; format 3,0; iter;; " ITERATIONS."; endif;
      else; rr[i] = exp(bm);
   endif;
endo;
if sw; rr = 1/rr[3 2 1]; endif;
if names[1]; /* print results: */
   if zed; "Point estimate undefined. ";;
      if cnv; if _midp; "Mid-p ";; else; "Fisher-p ";; endif;
         format /rdn 2,0; 100*_clevel;; "% "$+namerr[3-2*sw]$+":";;
         format /rds 10,4; rr[3-2*sw]'; else; ".";
      endif;
      else; "Mid-p rate ratio estimate and ";;
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