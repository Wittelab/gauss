/* This procedure does ordinal ML logistic regression using the
   cumulative-odds model.
 Inputs: x = design matrix,
         y = matrix of index outcome counts (do not include reference level),
         n = vector of totals,
         rep = vector of repetion counts (0 if none),
         offset = scalar or vector of offsets (0 if no offsets);
         bnames = vector of coefficient names
                  (set to 0 if no printed output is desired),
         ynames = vector of outcome names (may be 0 if no printed output).
Outputs: b = intercept estimates (one per index category)
             and slope coefficient estimates,
         bcov = inverse-information covariance matrix for vec(b),
         dev = residual deviance for model,
         df = residual degrees of freedom.
Globals (may be set by user from calling program; otherwise defaults are used):
      _maxit = maximum number of iterations (default is 30)
      _bcriter = convergence criterion for b (default is .001)
      _dcriter = convergence criterion for deviance (default is .1)
*/
proc (4) = lrordco(x,y,n,rep,offset,bnames,yname);
declare matrix _maxit = 30;
declare matrix _bcriter = .001; declare matrix _dcriter = .1;
local t0,est,nr,const,nc,nc1,nd,nd1,np,df,y0,yc,yv,xc,
      llsat,b,bcov,bcovr,bold,dev,devold,ind,iter,cnv,
      eta,p,mu,fit,i,w,xw,score;
t0 = date; bnames = 0$+bnames;
if sumc(rep); y = y.*rep; n = n.*rep; endif;
nr = rows(x); const = ones(nr,1);
/* No. covariates: */ nc = cols(x);
nc1 = nc + 1;
/* No. index outcome categories: */ nd = cols(y);
nd1 = nd + 1;
/* No. of parameters: */ np = nc + nd;
/* Degrees of freedom: */ df = nr*nd - np - sumc(n .eq 0);
if bnames[1]; print; format /rdn 3,0;
   "Proc lrordco.g: Ordinal logistic regression - cumulative-odds model";
   "     with ";; nd;;" index outcomes, and";; nc;;" regressors.";
   if rows(bnames) ne cols(x);
      "INPUT ERROR: No. names supplied not equal to no. of regressors."; end;
   endif;
endif;
/* Counts at reference level of y: */ y0 = n - sumr(y);
/* Saturated loglikelihood - sumc(ln(combin(n,y))): */
   llsat = sumall(y.*ln(y + (y .eq 0))) + y0'ln(y0 + (y0 .eq 0))
           /* - n'ln(n + n .eq 0) */ ;
/* Cumulative upper sums of y: */ yc = cumsumc(rev(y'))';
/* index matrix: */ ind = cumsumc(ones(nr,nd));
/* Paste an identity matrix alongside nd copies of each row of x
   and stack them up to make xc: */
   xc = (const.*.eye(nd))~(x[vecr(ind),.]);
local bp,binit,rss,inf,bsum,infsum;
/* Initialize b by taking an information-weighted average of the b from nd
   logistic regressions, one for each possible cutpoint for y: */
   bsum = 0; infsum = 0; bp = zeros(nc1,nd); i = 0;
   binit = logit(sumc(yc)./sumc(n))'|zeros(nc,nd);
   do until i ge nd; i = i + 1;
      { b,bcov,cnv } = lrquick(const~x,yc[.,i],n,offset,binit[.,i]);
      if cnv le 0; break; endif;
      bp[.,i] = b;
      inf = invpd(bcov[2:nc1,2:nc1]);
      infsum = infsum + inf;
      bsum = bsum + inf*bp[2:nc1,i];
   endo;
   if cnv le 0; /* Try simple initialization: */
      b = logit(sumc(yc)./sumc(n))|zeros(nc,1);
      else; b = bp[1,.]'|(invpd(infsum)*bsum);
   endif;
/* Reversed table of outcome counts: */ yv = (rev(y')')~y0;
/* Now fit cumulative-odds model by doing correlated logistic regression
   on the cumulative case-counts yc: */
   /* Initialize deviance, counter, convergence criterion: */
   dev = 0; devold = 1; iter = 0; cnv = 0;
   do until cnv or (iter ge _maxit); iter = iter + 1;
      /* linear predictors: */ eta = b[1:nd]' + (offset + x*b[nd1:np]);
      /* fitted cumulative probs: */ p = exp(eta); p = p./(1 + p);
      /* fitted cumulative case counts: */ mu = n.*p;
      /* fitted expectations for yv: */ fit = (mu~n) - (zeros(nr,1)~mu);
      /* new deviance: */ devold = dev; dev = 2*(llsat - sumall(yv.*ln(fit)));
      if bnames[1]; "Deviance at iteration ";;
         format /rds 4,0; iter;; ":";; format 12,3; dev;
      endif;
      if iter gt 1 and devold lt dev; dev = devold;
         /* step halve: */ b = (bold + b)/2; continue;
      endif;
      /* for each row of x,
         1) construct the covariance matrices for the cumulative counts,
            v = xpnd(vech(mu[ind[i,.],.]')) - mu[i,.]'p[i,.])
         2) construct the weight matrix w = invpd(1/v);
         3) xw = block i of xc times w: */
      xw = {}; i = 0;
      do until i ge nr; i = i + 1;
          w = invpd(1/(xpnd(vech(mu[ind[i,.],.]')) - mu[i,.]'p[i,.]));
          xw = xw~(w|(x[ind[i,.],.]'w));
      endo;
      score = xw*vecr((yc - mu)./(mu.*(1-p)));
      /* inverse information: */ bcov = invpd(xw*xc);
      /* Newton step: */ bold = b; b = b + bcov*score;
      cnv = prodc(prodc(abs(b - bold) .lt _bcriter))
            and (abs(devold - dev) lt _dcriter);
   endo;
if iter ge _maxit;
   "WARNING: No convergence after ";; format 4,0; iter;; " iterations";
endif;
if bnames[1];
    /* create names for intercepts: */
    bnames = (0$+"y = "$+ftocv(seqa(nd,-1,nd),2,0))|bnames;
    rreport(b,sqrt(diag(bcov)),bnames,yname,dev,-1,df,1);
    if (meanc(n) lt 5) and (df gt 0);
      "The mean denominator value (MDV) is ";; meanc(n);
      "The above deviance test of fit may not be valid if MDV < 5.";
    endif;
endif;
"Total run time: ";; etstr(ethsec(t0,date));
format 10,4;
retp(b,bcov,dev,df);
endp;
