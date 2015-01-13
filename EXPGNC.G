/* Gauss-Newton fitting of loglinear models to coarsened tables
   Inputs:
   x = complete design matrix
   aa = observed coarsened counts
   cn = numerator coarsening matrix 
        -- maps latent counts to observed counts aa
   cg = multinomial grouping matrix for coarsened (observed) counts
        -- maps observed counts to fixed margins (0 for Poisson)
   offset (0 if none)
   irand = indices in x of variables with random marginals (0 for Poisson)
   bp = prior mean for beta (0 for no prior)
   invp = prior inverse variance for beta (0 for no prior)
   bold = initial value for b
   Outputs:
   b = MLE or MPE of beta
   bcov = true inverse information for b
   eal = MLEs of complete counts under MCAR given sufficient stats for model
   ac = MLEs of complete counts given observed counts aa under model MCAR 
   lnpost = ln likelihood or ln posterior for b
   Globals: 
   _eminit = number of initial _em cycles (default = 200; used only if bold = 0)
   _gmaxit = maximum number of iterations (default = 100)
   _gmcrit = b convergence criterion
   _print = 1 for print, 0 for none
*/
proc (5) = expgnc(x,aa,cn,cg,offset,irand,bp,invp,bold);
declare _eminit = 200; declare _gmaxit = 100; declare _gcrit = .0001; 
declare _print = 1;
local t0,nb,nr,iter,mu,ac,eal,lnpost,score,inf,b,bcov,nit;
t0 = date;
nb = cols(x); nr = rows(x);
if _print; print; 
   "Gauss-Newton for loglinear modeling of coarsened tables";
   "cols x, rows input b = ";; format 1,0; nb;; rows(bold); format 8,3;
endif;
if rows(bold) ne nb;
   if _print; "Self-initialized with "$+ftocv(_eminit,1,0)$+" E-1-step cycles -- "; endif;
   bold = ln((x'(cn'(aa./sumr(cn))))./(sumc(x)+(sumc(x) .eq 0))); 
   iter = 0;
   do until iter ge _eminit; iter = iter+1;
      /* E step: */  mu = exp(offset+x*bold); 
                     ac = (cn'(aa./(cn*mu))).*mu;
      /* 1 step to M: */ 
         trap 1; bcov = invpd(moment(sqrt(mu).*x,0)+invp);
      if scalerr(bcov); 
        "EM INFORMATION NOT PD IN EXPGNC.G AT iter = "$+ftocv(iter,1,0);
         trap 0; end; 
      endif;
      bold = bold + bcov*(x'(ac-mu)-invp*(bold-bp)); 
   endo;
endif;
if irand[1]; 
   if rows(bp) eq nb; bp = bp[irand]; invp = invp[irand,irand]; endif;
   x = x[.,irand]; bold = bold[irand]; nb = cols(x); 
endif;  
iter=0; 
do until iter ge _gmaxit; iter = iter+1; 
   if irand[1]; {lnpost,eal,ac,score,inf} = postmc(aa,x,cn,cg,bold,offset,bp,invp); 
      else; {lnpost,eal,ac,score,inf} = postpc(aa,x,cn,bold,offset,bp,invp);
   endif;
   trap 1; bcov = invpd(inf);
   if scalerr(bcov); 
   "GN INFORMATION NOT PD IN EXPGNC.G AT iter = "$+ftocv(iter,1,0);
         trap 0; end; 
   endif;
   b = bold + bcov*score;
   if prodc(abs(b-bold) .lt _gcrit); break; endif; bold = b;
endo;
bold = b;
if iter ge _gmaxit;
   " EXPGNC WARNING: No convergence after iteration "$+ftocv(iter,1,0);
   elseif _print; ""$+ftocv(iter,1,0);;" GN iterations,";
endif;
if _print;
   " returning covariances from ";;
   if irand[1]; "multinomial";; else; "Poisson";; endif; 
   " information with ";;
   if rows(bp) eq nb; "prior."; else; "no prior."; endif;
   format /rds 8,3;
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(b,bcov,eal,ac,lnpost);
endp;
