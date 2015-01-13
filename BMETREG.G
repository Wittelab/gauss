/* This proc does Metropolis-sampling Bayesian risk or rate regression,
   with user-specified likelihood and prior.
   Objective is to estimate posterior mean, variance, and percentiles
   of the marginal posterior distributions for beta.
  INPUTS:
   x = covariate (design) matrix with np columns
   y = no. cases at each covariate level OR unknown-covariate indicator
   n = denominators for y at each covariate level OR observed margins
   nmc = no. of Markov chains
   crun = no. of discarded run-in cycles
   cuse = no. of cycles used after burn-in
   mets = no. of Metropolis steps
   bnames = names of x columns
   yname = name of y
   b = initial values for beta
   bcov = initial covariance for beta
   bp = prior location for beta
   t2 = squared prior scale
  Selected Intermediates:
   mcu = nmc*cuse = number of used cycles
   cl = crun+cuse = length of each chain
   yx = either y'x or selected rows of x
   b1, b2 = trial draws at second stage
   b = final draw at second stage, chosen from b1 & b2 with probabilities
       proportional to estimated posterior probablity density at b1 & b2
  OUTPUTS:
   betavg = mean of b over the last cuse trials of all chains
   betavgv = covariance matrix of b over the last cuse trials of all chains
   bcents = np by 7 matrix of 2.5, 5, 10, 50, 90, 95 & 97.5 percentiles of b
   bmat = nmc*cuse by np matrix of simulated b. THIS MAY BE A HUGE MATRIX!
  GLOBALS: _accept = scale factor for acceptance density (default is 1)
            proposal density is rescaled version of normal(b,bcov), see below
           _fobs = indicates fused-observation analysis (default 0)
            if _fobs, y distinguishes x-rows and splits x into xy and x and
                      n is the fused count, so sumc(y) = rows(x)/2 = rows(n)
*/
proc (4) = bmetreg(x,y,n,nmc,crun,cuse,mets,bnames,yname,b,bcov,bp,t2);
declare matrix _accept = 1;
declare matrix _fobs = 0;
local t0,np,yx,mcu,cl,invt2,cholbcov,binit,bmeans,bvars,bmat,bsq,nacc;
t0 = date; bnames = 0$+bnames; yname = 0$+yname;
if bnames[1]; "MCMC-METROPOLIS SAMPLING FROM EXACT POSTERIOR.";
   " - No. chains, run-in & used cycles, Metropolis steps:";;
   format /rds 1,0; nmc;; crun;; cuse;; mets; print; format 9,3;
endif;
/* Derived program constants: */
   /* total used samples: */ mcu = nmc*cuse;
   if mcu%40;
      "WARNING: NO. CHAINS TIMES NO. CYCLES NOT DIVISIBLE BY 40:";
      "PERCENTILES WILL BE INACCURATE.";
   endif;
   /* total chain length: */ cl = crun + cuse;
   np = cols(x);
/* For fused observations: */
   if _fobs; yx = selif(x,y); x = selif(x,1-y); else; yx = y'x; endif;
/* invert t2: */
   if cols(t2) eq 1; invt2 = eye(np)./t2;
      elseif cols(t2) eq np; invt2 = invpd(t2);
      else; "INPUT ERROR: t2 has "; dim(t2); end;
   endif;
/* Initialize for loop: */
   cholbcov = chol(bcov); binit = b + cholbcov*rndn(np,nmc);
   /* Chain-specific mean betas: */ bmeans = zeros(nmc,np);
   /* Chain-specific var betas: */ bvars = zeros(nmc,np);
   /* Accumulator for MCMC variance estimates: */ bsq = zeros(np,np);
   /* Matrix of generated betas: */ bmat = zeros(mcu,np);
   /* proposal acceptance count: */ nacc = 0;
local j,rnb,mjump,u,i,itot,mi,b1,b2,lpr;
/* Outer (j) loop goes across Markov chains: */
"GENERATING CHAIN ";; j = 0;
do until j ge nmc; j = j + 1; format 1,0; j;; format 9,3;
   /* Random nos. transformed to proposal deviations (Gelman et al. p335): */
      mjump = _accept*((2.4/sqrt(np))*cholbcov)'rndn(np,cl*mets);
   /* Uniform variates for choice between b1 and b2: */ u = ln(rndu(cl,mets));
   /* Inner (i) loop through one chain: */ b = binit[.,j]; i = 0;
   do until i ge cl; i = i + 1;
      /* Metropolis cycle: */ b1 = b; mi = 0;
         do until mi ge mets; mi = mi + 1;
            /* normal proposal density with prior variance rescaled: */
               b2 = b1 + mjump[.,(mi-1)*cl+i];
            /* Compute ln posterior ratio, b2 vs. b1: */
               lpr = lnpratio(yx,x,y,n,b2,b1,bp,invt2);
            /* check u<lpr to select b1 vs. b2: */
            /* Metropolis step -- take b2 if post prob(b2) > post prob(b1);
               if not, take b2 with probability ppr: */
            if u[i,mi] lt lpr; b1 = b2; nacc = nacc + 1; endif;
         endo;
         b = b1;
    /* After run-in, begin accumulating results: */
       if i gt crun;
          /* index of this entry: */ itot = (cuse*(j-1) + (i-crun));
          bmat[itot,.] = b';
          bmeans[j,.] = bmeans[j,.] + b';
          bvars[j,.] = bvars[j,.] + (b.*b)';
          bsq = bsq + b*b';
       endif;
   endo;
endo;

print; print; format 10,5; print;
local betavg,betavgv,bmeancov;
bmeans = bmeans/cuse;
betavg = meanc(bmeans);
bmeancov = (bmeans - betavg')'(bmeans - betavg');
betavgv = bsq/mcu - betavg*betavg' + bmeancov/nmc;

if bnames[1];
   "MCMC results with prior variances = "; t2';
   format 9,3; "Acceptance rate:";; nacc/(nmc*cl*mets);

   /* slow & inaccurate mode finding:
      local bmode,h,bins,nbins,cbin;
      i = 0; bmode = zeros(np,1); h = .05;
      bins = seqa(-3,h,1 + 8/h); nbins = rows(bins);
      do until i ge np; i = i + 1; j = 0; cbin = zeros(nbins,1);
         /* looped to avoid memory overflow: */
         do until j ge nbins; j = j + 1;
          cbin[j] = sumc((bmat[.,i] .gt bins[j]).*(bmat[.,i] .le (h+bins[j])));
         endo;
         bmode[i] = bins[maxindc(cbin)];
      endo;
      "Mode estimate and antilog:";; bmode;; exp(bmode);
   */

   /* "SDs of sim. betas computed two ways:";
      sqrt(diag(betavgv))~sqrt((mcu-1)/mcu)*stdc(bmat); */

   /* Compute Gelman-Rubin convergence criterion (Stat Sci 1992): */
    local m1m,c1c,bm2,w,bn,v,covbvbm2,covbvbm,vv,r;
    if nmc gt 1;
       m1m = (nmc+1)/nmc; c1c = (cuse-1)/cuse; bm2 = bmeans.*bmeans;
       bvars = (bvars - cuse*bm2)/(cuse-1); w = meanc(bvars);
       bn = nmc*(meanc(bm2) - betavg.*betavg)/(nmc-1);
       v = c1c*w + m1m*bn;
       covbvbm2 = (meanc(bvars.*bm2) - w.*meanc(bm2))/(nmc-1);
       covbvbm = (meanc(bvars.*bmeans) - w.*betavg)/(nmc-1);
       vv = ((c1c*stdc(bvars))^2)/nmc + 2*((m1m*bn)^2)/(nmc-1)
            + 2*m1m*c1c*(covbvbm2 - 2*betavg.*covbvbm);
       r = v./sqrt(w.*(v - vv./v));
       else; r = zeros(np,1);
    endif;

   local border,icents,bcents,wald,si,mask8,fmt8;
   let icents = .025 .05 .1 .5 .9 .95 .975; icents = icents*mcu;
   bcents = zeros(np,7); i = 0;
   do until i ge np; i = i + 1;
      border = sortc(bmat[.,i],1);
      bcents[i,.] = border[icents]';
   endo;

   let mask8[1,8] = 0 1 1 1 1 1 1 1;
   let fmt8[8,3] = "-s" 8 8  "lf," 8 3  "lf," 8 3  "lf," 8 3
                    "lf," 8 2  "lf," 8 2  "lf," 8 2  "lf" 8 2;
   "Gelman-Rubin criteria (0 if not computed because only one chain),";
   " prior means, means of simulated betas, medians of exp(betas),";
   " and geomean and Wald limits for simulated exp(betas):";
   wald = exp(betavg+(0~(-1.96)~1.96).*sqrt(diag(betavgv)));
   call printfm(bnames~r~bp~betavg~exp(bcents[.,4])~wald,mask8,fmt8);
   print;
   "Serial correlation of estimates: "; format /rds 5,2;
   si = seqa(1,1,rows(bmat)-1); i = 0;
   do until i ge np; i = i + 1;
      corrp(bmat[si,i],bmat[si+1,i]);;
   endo; print; print;
    let mask8[1,8] = 0 1 1 1 1 1 1 1;
   /*
    let fmt8[8,3] = "-s" 8 8  "lf," 8 3  "lf," 8 3  "lf," 8 3
                    "lf," 8 3  "lf," 8 3  "lf," 8 3  "lf" 8 3;
    " 2.5, 5 , 10, 50, 90, 95, & 97.5 % of simulated betas:";
    call printfm(bnames~bcents,mask8,fmt8); print;
   */
    let fmt8[8,3] = "-s" 8 8  "lf," 8 2  "lf," 8 2  "lf," 8,2
                    "lf," 8 2  "lf," 8 2  "lf," 8 2  "lf" 8 2;
    " 2.5, 5 , 10, 50, 90, 95, & 97.5 % of simulated exp(betas):";
    call printfm(bnames~exp(bcents),mask8,fmt8); print;
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(betavg,betavgv,bcents,bmat);
endp;

proc lnpratio(yx,x,y,n,b2,b1,bp,invt2);
/* fused-exponential loglikelihood, using y to indicate the unknown covariate,
   yx and x as the submatrices of the original x with and without y=1,
   and n as the observed margin: */
   /* To reject crossovers:
      if sign(b2[7]+b2[8]) ne sign(b2[7]); retp(-1000); endif; */
   local mu2,mu1; mu2 = exp(yx*b2)+exp(x*b2); mu1 = exp(yx*b1)+exp(x*b1);
retp(
/* logistic llr: yx*(b2-b1) - n'ln((1+exp(x*b2))./(1+exp(x*b1))) */
/* exponential llr: yx*(b2-b1) - n'(exp(x*b2) - exp(x*b1)) */
/* fused-exponential llr: */ n'ln(mu2./mu1) - sumc(mu2 - mu1)
/* normal ln prior ratio: */
-((b2-bp)'invt2*(b2-bp) - (b1-bp)'invt2*(b1-bp))/2
);
endp;
