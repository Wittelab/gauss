/* This procedure does logistic regression via the Flanders-Greenland
   pseudo-likelihood (Stat Med 1991) for nested and two-stage studies.
 Inputs: x = design matrix,
         y = vector of case counts,
         n = vector of totals (if scalar, all totals will be set to this n),
         s = vector giving the sampling-stratum variable for each record,
         u = pseudolikelihood weight (inverse of sampling probability)
             for each record (it is assumed u is constant within strata),
         rep = vector of repetion counts (0 if none),
         const = 0 if no intercept (constant),
         offset = vector of offsets (0 if no offsets),
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
         yname = scalar name of outcome (may be 0 if no printed output).
 Outputs: b = coefficient estimates,
          bcovr = robust ("sandwich") covariance matrix for b,
          dev = residual pseudo-deviance for model,
          df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _binit = vector of initial values for b (default is 0; if user-specified
               and constant is requested, it must include a value for the
               constant as its first element)
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
*/
proc (4) = lrpseudo(x,y,n,s,u,rep,const,offset,bnames,yname);
declare matrix _binit = 0; declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .1;
local t0,su,si,uu,np,nr,yu,nu,neq0,nneq0,df,llsat,b,bold,dev,devold,
      iter,cnv,eta,p,mu,w,infi,res,scores,bcovr;
t0 = date;
su = unique(s,1);
/* stratum indicators: */ si = s .eq su';
/* mean stratum weights: */ uu = (si'u)./sumc(si);
if sumc(sumr(uu' .eq u.*si) .ne 1);
   "ERROR IN LRPSEUDO.G: Weights (u) must be equal within strata (s).";
end; endif;
if bnames[1]; print;
   " Proc lrpseudo.g: Pseudo-ML logistic regression for data with";
   "                  inverse-sampling fraction weights."; print;
   if rows(bnames) ne cols(x);
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
endif;
nr = rows(x);
if sumc(rep) ne 0; y = y.*rep; n = n.*rep; endif;
/* Include constant if requested: */ if const; x = ones(nr,1)~x; endif;
/* No. parameters: */ np = cols(x);
/* Flag empty covariate patterns: */ neq0 = n .eq 0; nneq0 = n + neq0;
/* Degrees of freedom: */ df = nr - np - sumc(neq0);
yu = y.*u; nu = n.*u;
/* Saturated loglikelihood + n'ln(nneq0) - sumc(ln(combin(n,y))): */
   llsat = yu'ln(y + (y.==0)) + (nu - yu)'ln(n - y + (y.==n));
/* Initialize beta, deviance, counter,convergence indicator: */
   if rows(_binit) eq np; b = _binit; else; b = zeros(np,1); endif;
   dev = 0; devold = 1; iter = 0; cnv = 0;
/* Begin iterative reweighted LS estimation */
do until cnv or (iter ge _maxit); iter = iter + 1;
   /* linear predictors: */ eta = offset + x*b;
   /* expected proportions: */ p = 1./(1 + exp(-eta));
   /* expected counts: */ mu = n.*p;
   devold = dev; dev = 2*(llsat - (yu'eta + nu'ln(nneq0 - mu)));
   if bnames[1]; "Pseudo-deviance at iteration ";;
      format /rds 4,0; iter;; ":";; format 12,3; dev; endif;
   if iter gt 1 and devold lt dev; dev = devold;
      /* step halve: */ b = (bold + b)/2; continue;
   endif;
   /* u-weight vector: */ w = u.*mu.*(1-p);
   /* covariance of b: */ infi = invpd((w.*x)'x);
   /* u-residual: */ res = u.*(y - mu);
   /* Newton step: */ bold = b; b = b + infi*(x'res);
   cnv = prodc(abs(b - bold) .lt _bcriter) and (abs(devold - dev) lt _dcriter);
endo;
if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
/* matrix of stratum-specific score sums: */ scores = si'(x.*res);
/* bcovr successively built into robust covariance matrix: */
      bcovr = (x'((u.*u.*(y + p.*(mu - 2*y))).*x));
   /* 1-1/uu is the finite-population correction factor
      and si'n is the vector of stratum counts; 
      if si = identity, the following reduces to (x.*res./u)'(x.*res): */
      bcovr = bcovr - scores'(((1-1/uu)./(si'nneq0)).*scores);
   bcovr = infi*bcovr*infi;
if bnames[1]; print; "With specification-robust variance:";
   rreport(b,sqrt(diag(bcovr)),bnames,yname,-1,-1,df,1);
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(b,bcovr,dev,df);
endp;
