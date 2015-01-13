/* This procedure does kernel-smoothed logistic regression via local likelihood
 Inputs: xf = matrix of points to be fitted -- set to 0 if xf = x,
         is = indices of columns of xf that define kernels
              (if 0, it will be assumed all should be used)
         h = rows(xf) by rows(is) matrix of half-bandwidths for elements
             of xf[.,is] (if vector conformable with xf[.,is] or scalar,
              it will be expanded to matrix),
         x = design matrix,
         y = vector of case counts,
         n = vector of totals (if scalar, all totals will be set to this n),
         rep = vector of repetion counts (0 if none),
         const = 1 if you wish the program to add a constant to x, 0 if not;
                 if the first j (>1) columns of x represent a complete set
                 of indicators for one covariate, const = j will result in
                 weighted-mean centering of the remaining cols(x)-j covariates,
                 where the weights are those obtained from the final iteration,
                 which "approximately diagonalizes" the covariance submatrix
                 of the j indicator-coefficient (intercept) estimates.
         offset = vector of offsets (0 if no offsets),
         fnames = vector of names for fitted points
                  (if 0, the names "Pt 1", "Pt 2", etc. will be used),
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
         yname = scalar name of outcome (may be 0 if no printed output).
 Outputs: bl = matrix of local coefficient estimates, one column
               for each row of xf (with intercept first if const = 1)
          bcovr = vertical concatenation of the robust ("sandwich") covariance
                  matrices for the columns of b,
          lf = fitted logits corresponding to xf,
          lcov = covariance matrix for lf,
          The following are computed only if xf = x (otherwise 0 returned):
          dev = residual deviance for model,
          rss = residual weighted sum of squares,
          df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _clevel = confidence level expressed as a proportion (default is .95)
      _fit = 1 for iteration reporting and tests of fit, 0 if not (default is
             1, but no tests provided if mean of case or noncase count is <3)
      _mis = 0 if no missing values, 1 if complete-case analysis
      _mcode = missing-value scalar or column vector with element
               for each column of x (default -99999; ignored if _mis = 0)
      _binit = vector of initial values for b (default is weighted-least
               squares estimates; if user-specified and constant is
               requested, it must include a value for the constant as its
               first element)
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
      _specrob = 1 if you want an additional version of output
                 in which the Huber specification-robust "sandwich"
                 covariance matrix estimate (bcovr) is used instead of
                 the usual inverse-information estimate (default is 0)
      _bias = 1 if you want an additional printout of odds ratios
              computed with the Cox-Snell ML bias correction (default is 0)
*/
proc (7) = logregs(xf,is,h,x,y,n,rep,const,offset,fnames,bnames,yname);
declare matrix _mis = 0; declare matrix _mcode = -99999;
declare matrix _binit = 0; declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .1;
declare matrix _fit = 1; declare matrix _clevel = .95;
local t0,neq0,nr,mcc,np,llsat,b,bold,bcov,bcovr,dev,devold,df,iter,cnv,
      eta,p,mu,w,rss;
t0 = date; bnames = 0$+bnames; yname = 0$+yname;
if sumc(rep); y = y.*rep; n = n.*rep; endif;
/* Delete empty covariate patterns: */ neq0 = n .eq 0;
   if _mis eq 1; /* also delete records with missing regressor values: */
      neq0 = neq0 .or sumr(x .eq _mcode'); endif;
   if sumc(neq0); x = delif(x,neq0); y = delif(y,neq0);
      if rows(n) eq rows(neq0); n = delif(n,neq0); endif;
      if rows(offset) eq rows(neq0); offset = delif(offset,neq0); endif;
   endif;
nr = rows(x);
if rows(n) eq 1; n = n*ones(nr,1); endif;
if bnames[1]; print;
   "PROC LOGREGS.G: Kernel-smoothed logistic regression for "$+yname$+".";
   if rows(bnames) ne cols(x);
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
   mcc = meanc(y)|meanc(n-y);
   print; "Nos. of cases and noncases :";; format /rds 12,0; nr*mcc';
endif;
/* Add constant if requested: */ if const eq 1; x = ones(nr,1)~x; endif;
/* No. parameters: */ np = cols(x);
if rank(x) lt np; "DESIGN MATRIX RANK DEFICIENT IN LOGREGS.G";
   retp(0,0,0,0,0,0,0); endif;
/* Restricted degrees of freedom: */ df = nr - np;
/* Saturated loglikelihood + n'ln(n) - sumc(ln(combin(n,y))): */
   llsat = y'ln(y + (y.==0)) + (n - y)'ln(n - y + (y.==n));

/* Begin local logistic fitting: */
   local gfit,nfit,k,bl,pu,bcovi,u,u0,yl,nl,lf,lcov,devl,rssl,llsatl,vk,j,vj;
gfit = rows(xf) eq 1 and cols(xf) eq 1 and xf[1,1] eq 0; nfit = rows(xf);
if not fnames[1]; fnames = 0$+"Pt "$+ftocv(seqa(1,1,nfit),1,0); endif;
if rows(fnames) ne nfit;
   "ERROR IN LOGREGS.G: No. fit names not equal to no. fit points"; end; endif;
if is[1] eq 0; is = seqa(1,1,np-const); endif;
if const; is = 1 + is; endif;
if gfit; xf = x; elseif const; xf = ones(rows(xf),1)~xf; endif;
if rows(h) eq 1; h = h[ones(nfit,1),.]; endif;
bl = zeros(np,nfit); pu = zeros(rows(x),nfit);
bcovi = zeros(np*nfit,np); bcovr = bcovi;
lf = zeros(nfit,1); lse = lf; k = 0;
DO until k ge nfit; k = k + 1;
/* generalized Epanechnikov kernel weights = positive part of
   no. parameters minus squared Euclidean distances
   from fitted point in h-units: */
/* selif(x[.,2],abs(xf[k,2]-x[.,2]) .lt h[k]); */
   u = max(rows(is) - sumr(((xf[k,is]-x[.,is])./h[k,.])^2),0); u0 = u .eq 0;
   if sumc(1-u0) le np;
      "TOO FEW OBSERVATIONS IN KERNEL OF "$+fnames[k]; continue; endif;
   yl = u.*y; nl = u.*n;
   if _fit*bnames[1];
      "Kernel-weighted means for cases and noncases at "$+fnames[k];;":";
      format 9,3; $bnames[1+const:np]';
      meancw(x[.,1+const:np],y.*u)';
      meancw(x[.,1+const:np],(1-y).*u)';
   endif;
/* Saturated local loglikelihood + n'ln(n) - sumc(ln(combin(n,y))): */
   llsatl = yl'ln(yl + (yl.==0)) + (nl - yl)'ln(nl - yl + (yl.==nl));
/* Initialize beta, deviance, convergence indicator, counter: */
   if rows(_binit) eq np; b = _binit;
      else; /* start with WLS estimates: */ format 9,3;
      p = (yl + sumc(yl)/sumc(nl) + .01)./(nl + 1.02);
      w = nl.*p.*(1-p);
      trap 1; bcov = invpd(moment(sqrt(w).*x,0));
      if scalerr(bcov); "SINGULARITY INITIALIZING LOGREGS.G AT "$+fnames[k];
         continue; endif; trap 0;
      b = bcov*(x'(w.*(logit(p)-offset)));
   endif;
dev = 0; cnv = 0; iter = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictors: */ eta = offset + x*b;
   /* expected proportions: */ p = 1./(1 + exp(-eta));
   /* expected counts: */ mu = nl.*p;
   if sumc(nl-mu+u0 .eq 0);
      "DIVERGENCE IN LOGREGS.G AT "$+fnames[k]; u[1] = -1; break; endif;
   devold = dev; dev = 2*(llsat - (yl'eta + nl'ln(nl - mu + u0)));
   if _fit*bnames[1]; "Deviance at iteration ";;
      format /rds 4,0; iter;; ":";; format 12,3; dev; endif;
   if iter gt 1 and devold lt dev; dev = devold;
      /* step halve: */ b = (bold + b)/2; continue;
   endif;
   /* weight vector: */ w = mu.*(1-p);
   /* covariance of b: */ trap 1; bcov = invpd(moment(sqrt(w).*x,0));
      if scalerr(bcov); "SINGULARITY IN LOGREGS.G AT "$+fnames[k];
         u[1] = -1; break; endif; trap 0;
   /* Newton step: */ bold = b; b = b + bcov*(x'(yl-mu));
   cnv = prodc(abs(b - bold) .lt _bcriter) and (abs(devold - dev) lt _dcriter);
endo;
if u[1] lt 0; continue; endif;
/* weighted sum of squared residuals: */ rss = (yl-mu)'((yl-mu)./(w+u0));
/* pseudo-residual: e = (yl - mu)./w; */
/* pseudo-observed logit: z = eta + e; */
/* note that b = bcov*wx'z; */
if const gt 1 and np gt const; local kc; kc = seqa(const+1,1,np-const);
   /* center all covariates after the first j, then recompute bcov & b: */
   x[.,kc] = x[.,kc] - (w/sumc(w))'x[.,kc];
   bcov = invpd(moment(sqrt(w).*x,0)); b = bcov*(x'(w.*eta + yl - mu));
endif;
/* local fitted logit: */ lf[k] = xf[k,.]*b;
   bl[.,k] = b; pu[.,k] = p.*u;
   vk = 1+seqa((k-1)*np,1,np); bcovi[vk,.] = bcov;
if iter ge _maxit; "WARNING: No convergence after iteration";;
   format 1,0; iter;; " at point "$+fnames[k];
endif;
if _fit*bnames[1]; print; "Using specification-robust variances:";
   rreport(b,sqrt(diag(bcov)),bnames,yname,dev,
           rss*(-1)^(1-_fit),df,1+_fit*(minc(mcc) ge 3)); endif;
ENDO;
/* compute the covariance of lf: */
   lcov = zeros(nfit,nfit); j = 0;
   do until j ge nfit; j = j + 1; vj = 1+seqa((j-1)*np,1,np);
      bcovr[vj,.] = (bcovi[vj,.]*x'(((1-pu[.,j]).*y.*(1-pu[.,j])
                    + pu[.,j].*(n-y).*pu[.,j]).*x)*bcovi[vj,.]);
      lcov[j,j] = xf[j,.]*bcovr[vj,.]*xf[j,.]';
      k = j;
      do until k ge nfit; k = k + 1; vk = 1+seqa((k-1)*np,1,np);
         lcov[j,k] = (xf[j,.]*bcovi[vj,.]*x'(((1-pu[.,j]).*y.*(1-pu[.,k]) +
                      pu[.,j].*(n-y).*pu[.,k]).*x)*bcovi[vk,.]*xf[k,.]');
         lcov[k,j] = lcov[j,k];
      endo;
   endo;
if gfit; mu = n./(1 + exp(-lf));
   dev = 2*(llsat - (y'lf + n'ln(n - mu)));
   rss = (y-mu)'((y./mu-1)./(1-mu./n));
   else; dev = 0; rss = 0; df = 0;
endif;

if bnames[1]; "Using specification-robust variances:";
   if _fit*gfit;
      rreport(lf,lse,fnames,yname,dev,rss,df,1+(minc(mcc) ge 3));
      "The means of the case & noncase counts are ";; format 6,1; mcc';
      if minc(mcc) lt 5;
         "The minimum of the mean case & mean noncase counts is less than 5.";
         if minc(mcc) ge 3;
            "Hence the above deviance and rss tests of fit may be invalid.";
         endif;
      endif;
      else; rreport(lf,lse,fnames,yname,-1,-1,df,1);
   endif;
   if _fit; "Total run time: ";; etstr(ethsec(t0,date)); endif; print;
endif;
retp(bl,bcovr,lf,lcov,dev,rss,df);
endp;

