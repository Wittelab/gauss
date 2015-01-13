/* EM fitting of loglinear models to coarsened tables
   Inputs:
   x = complete-data design matrix
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
   bcove = true inverse information for b
   eal = MLEs of complete counts under MCAR given sufficient stats for model
   ac = MLEs of complete counts given observed counts aa under model MCAR 
   lnposte = ln likelihood or ln posterior for b
   Globals: 
   _emaxit = maximum number of full E-M cycles (default 2000)
   _emcrit = actual b convergence criterion
   _fmaxit = maximum number of iterations for each M step (default 1)
   _fmcrit = b convergence criterion for each M step
*/
proc (5) = expemc(x,aa,cn,cg,offset,irand,bp,invp,bold);
declare _emaxit = 2000; declare _emcrit = .000001;
declare _fmaxit = 1; declare _fmcrit = .001;
local t0,nb,nr,tit,iter,mu,ac,jter,boldm,b,bcove,nit;
t0 = date;
print; "EM allowing up to "$+ftocv(_fmaxit,1,0)$+" iterations per M step --";
nb = cols(x); nr = rows(x); 
"cols x, rows input b = ";; format 1,0; nb;; rows(bold); format 8,3;
if rows(bold) ne nb; "Self-Initializing.";
   bold = ln((x'(cn'(aa./sumr(cn))))./(sumc(x)+(sumc(x) .eq 0))); endif;
tit=0; iter=0;
do until iter ge _emaxit; iter = iter+1;
   /* E step: */  mu = exp(offset+x*bold); 
                  ac = (cn'(aa./(cn*mu))).*mu;
   /* M step: */ jter = 0; boldm = bold;
   do until jter ge _fmaxit; jter = jter+1; tit = tit+1;
         trap 1; bcove = invpd(moment(sqrt(mu).*x,0)+invp);
      if scalerr(bcove); 
        "EM INFORMATION NOT PD IN EXPEMC.G AT iter = "$+ftocv(iter,1,0);
         trap 0; end; 
      endif;
       b = boldm + bcove*(x'(ac-mu)-invp*(boldm-bp)); 
       if prodc(abs(b-boldm) .lt _fmcrit); break; endif; boldm = b;
   endo;
   if prodc(abs(b-bold) .lt _emcrit); break; endif; bold = b;
endo;
bold = b;
if iter ge _emaxit;
   " EXPEMC WARNING: No convergence after iteration "$+ftocv(iter,1,0);
   else;
   " convergence in "$+ftocv(iter,1,0)$+" EM cycles and "$+ftocv(tit,1,0)$+" inversions";
endif;
local lnposte,eal,score,inf; " returning covariances from ";;
if irand[1]; "multinomial information";;
   if rows(bp) eq nb; bp= bp[irand]; invp = invp[irand,irand]; endif;
   {lnposte,eal,ac,score,inf} = postmc(aa,x[.,irand],cn,cg,b[irand],offset,bp,invp);
   else; "Poisson information";;
   {lnposte,eal,ac,score,inf} = postpc(aa,x,cn,b,offset,bp,invp); 
endif;
if rows(bp) eq nb; " + prior."; else; ", no prior."; endif;
bcove = invpd(inf); format /rds 8,3;
"Total run time: ";; etstr(ethsec(t0,date));
retp(b,bcove,eal,ac,lnposte);
endp;
