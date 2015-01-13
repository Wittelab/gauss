load dat[32,6] = /stat/egret/data/williams.dat;
/* The 6 cols of dat are denom  outcome group grp1 grp2  matching
   EGRET random-effects normal-logistic regression yields
   constant = 2.624 (SE = .483), group b = -1.080 (SE = .626) for these data.
 */
x = dat[.,3]; y = dat[.,2]; n = dat[.,1]; id = dat[.,6];
let bnames = "group"; let yname = "n defect";
output file = temp.out reset;
nr = rows(x); format 10,3;
/* No repetition counts: */ a = 0;
/* No offset: */ offset = 0;
/* I want the constant in the regression: */ const = 1;
t = 0;
 { b,bcov,bcovr,dev,rss,df } =
               logreg(x,y,n,a,const,offset,bnames,yname);
{ b,bcovr,r,s,rss,df } = lrgee(x,y,n,x,t,const,offset,bnames,yname);
end;
/* This proc fits logistic regression models to clustered (correlated)
   binary-outcome data with possibly unequal-sized clusters, using either
   exchangeable or first-order autoregressive correlation structures, via
   a modification of the Liang-Zeger generalized estimating equation method.
 Inputs: x = design matrix,
         y = vector of case counts,
         n = vector of totals,
         id = cluster identification variable (e.g., individual id no. for
              longitudinal data analysis)
         t = vector of times of measurements (for autoregressive correlation)
             - may also be included in x
             - set to zero for exchangeable (constant) correlation structure,
         const = 0 if no intercept (constant),
         offset = vector of offsets (0 if no offsets),
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
         yname = scalar name of outcome (may be 0 if no printed output).
Outputs: b = coefficient estimates,
         bcovr = robust ("sandwich") covariance matrix for b,
         r = estimated correlation between adjacent cluster observations,
         s = Liang-Zeger estimate of estimate for scale factor (inverse of phi),
             same as ordinary residual sum of squares,
         rss = residual correlation-matrix-weighted sum of squares,
         df = residual degrees of freedom.
*/
proc (6) = lrgee(x,y,n,id,t,const,offset,bnames,yname);
/* Max. iterations, convergence criteria: */
declare matrix _bcriter = .001;
declare matrix _rcriter = .1;
declare matrix _maxit = 30;
local t0,np,nr,b,bold,r,s,rnew,rss,evar,df,neq0,is,cs,nc,rdf2,iter,cnv,
      eta,p,mu,w,res,i,ind,sw,pres,inf,score,swx,invr,presi,temp,bcovr,mdv;
/* local df,neq0,is,cs; */
t0 = date;
if bnames[1]; print;
   " Proc lrgee.g: logistic regression for clustered (correlated) outcomes.";
   " via modified Liang-Zeger generalized estimating equations.";
   " WARNING: CLUSTER MEMBERS MUST HAVE ADJACENT RECORDS - ";
   "  SORT ON CLUSTER ID BEFORE USE!"; print;
endif;
/* No. rows: */ nr = rows(x);
/* Include constant if requested: */ if const; x = ones(nr,1)~x; endif;
/* No. parameters: */ np = cols(x);
/* No. empty covariate patterns: */ neq0 = n .eq 0;
/* Degrees of freedom: */ df = nr - np - sumc(neq0);
/* Identify cluster starts and sizes: */ { is,cs } = changept(id);
/* No. clusters: */ nc = rows(is);   is~cs;
/* two times the df for correlation: */ rdf2 = cs'(cs-1) - 2*np;
/* Initialize beta, correlation, counter: */
   b = zeros(np,1); rnew = 0; iter = 0; cnv = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit);
   iter = iter + 1;
   /* save old values: */ bold = b; r = rnew;
   /* linear predictors: */ eta = offset + x*b;
   /* expected proportions: */ p = 1./(1 + exp(-eta));
   /* expected counts: */ mu = n.*p;
   /* residuals: */ res = y - mu;
   /* square-root of weight vector: */ sw = sqrt(mu.*(1-p));
   /* ordinary Pearson residuals: */ pres = res./(sw + neq0);
   /* scale-factor estimate: */ s = (pres'pres)/df;
   /* loop through clusters to build score and information, and update r: */
   i = 0; inf = 0; score = 0; evar = 0; rnew = 0; rss = 0;
   do until i ge nc; i = i + 1;
      /* indices of records for cluster i: */ ind = seqa(is[i],1,cs[i]);
      /* invert correlation matrix for cluster: */
         invr = invpd(s*diagrv(r*ones(cs[i],cs[i]),ones(cs[i],1)));
      swx = sw[ind].*x[ind,.];
/* dim(ind); dim(invr); */
      temp = invr*swx;
      presi = pres[ind];
      /* information contributions: */ inf = inf + swx'temp;
      /* score contributions: */ score = score + swx'(invr*presi);
      temp = presi'temp;
      /* empirical-variance contributions: */ evar = evar + temp'temp;
      /* correlation contributions: */ rnew = rnew + sumall(presi*presi');
      /* RSS: */ rss = rss + presi'(invr*presi);
   endo;
   rnew = (rnew - pres'pres)/(rdf2*s);
   /* Newton step: */ b = b - invpd(inf)*(x'res);
   if bnames[1]; "RSS & scale factor at iteration ";;
      format /rdn 4,0; iter;; ":";; format 11,3; rss;; s; endif;
   cnv = prodc(abs(b - bold) .lt _bcriter) and (abs(rnew - r) lt _rcriter);
endo;
/* invert information: */ inf = invpd(inf);
/* robust ("sandwich") covariance: */ bcovr = inf*evar*inf;
if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
if bnames[1];
   print; "With specification-robust variance:";
   rreport(b,sqrt(diag(bcovr)),bnames,yname,-1,rss,df,1);
   mdv = meanc(n);
   if (mdv lt 5) and (df gt 0);
      "The mean denominator value (MDV) is ";; mdv;; ";";
      "    the above rss test of fit may be invalid if this is <5.";
   endif;
"Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(b,bcovr,r,s,rss,df);
endp;
