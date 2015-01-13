/* ln-posterior evaluation with coarsened loglinear multinomial counts, 
   assuming coarsening takes place within but not across the multinomials.
   Inputs:
   aa = coarsened numerators (fused multinomial counts)
   x = complete-data design matrix
   m = denominators (multinomial totals)
   cn = numerator coarsening matrix 
        -- maps latent counts to observed counts
   cg = multinomial grouping matrix for coarsened (observed) counts
        -- maps observed counts to fixed margins
   b = parameter value at which to evaluate posterior
   bp, invp = normal-prior mean and covariance (0 for no prior) 
   Outputs: 
   ln posterior (ln likelihood if no prior)
   eal = expected complete (latent) data 
   ac = expected complete data given observed aa
   ln posterior derivative (score if no prior)
   expected negative ln posterior Hessian (expected information if no prior)
   Global: _sho = no printed output, 1 = show SEs, 
                  2 = show matrices and SEs and stop
*/
proc (5) = postmc(aa,x,cn,cg,b,offset,bp,invp);
declare matrix _sho = 0;
local m,cm,mu,cmu,mmu,p,eal,cp,eaa,pw,ac,vw,vc,infaa,dn,dm,lnlik,score,lnprior;
/* multinomial totals: */ m = cg*aa;
/* maps latent counts to fixed margins */ cm = cg*cn;
/* odds factors: */ mu = exp(offset+x*b); 
/* coarsened and marginal factors: */ cmu = cn*mu; mmu = cm*mu; 
/* latent probabilities and expected counts: */ p = (cm'(1/mmu)).*mu; eal = (cm'm).*p; 
/* coarsened probabilities and expectations: */ cp = cn*p; eaa = (cg'm).*cp;
/* probs and covariances of latent al given aa: */ 
   pw = (cn'(1/cmu)).*mu; ac = (cn'aa).*pw;
   vw = diagm(ac)-pw.*(cn'(aa.*cn)).*pw'; 
/* complete-data var */ vc = diagm(eal) - p.*(cm'(m.*cm)).*p'; 
/* actual beta inf given observed: */ infaa = x'((vc - vw)*x); 
/* each formula gives same ln likelihood: */ 
   if 0; lnlik = aa'ln(cmu) - m'ln(mmu); else; lnlik = aa'ln(cp); endif;
/* each formula gives same lnlik gradient: */ score = x'(ac - eal);
if _sho; local vaa,inf,lgrad; print;
   /* cov(aa): */ vaa = diagm(eaa)-cp.*(cg'(m.*cg)).*cp'; 
   if _sho gt 1; cp'; format 5,2;
      "cg, cg'cg:";; cg; print; cg'cg; print; "cov observed counts:";; vaa; print; 
      "cn, cn'cn:";; cn; print; cn'cn; print; "cov latent given observed:";; vw; print;
      "cm, cm'cm:";; cm; print; cm'cm; print; "cov latent (complete) counts:";; vc; print; 
      end;
   endif;
   /* derivatives of cmu, mmu, lnlik: */ dn = cn*(mu.*x); dm = cm*(mu.*x);
   format 8,3; lgrad = dn'(aa./cmu) - dm'(m./mmu); 
   "gradient:";; lgrad'; "score   :";; score';
   /* -hessian(lnlik): */ inf = dn'(dn./cmu) - dm'(dm./mmu);
   if rank(inf) ge rows(inf); "Naive SEs from -hess  :";; sqrt(diag(invpd(inf)))'; endif;
   inf = x'cn'(vaa*cn*x);
   if rank(inf) ge rows(inf); "Naive SEs from cov obs:";; sqrt(diag(invpd(inf)))'; endif;
   inf = x'(vc*x);
   if rank(inf) ge rows(inf); "complete-data SEs     :";; sqrt(diag(invpd(inf)))'; endif;
   if rank(infaa) ge rows(infaa); "observed-data SEs     :";; sqrt(diag(invpd(infaa)))'; endif; 
endif;
lnprior = -((b-bp)'invp*(b-bp))/2;
retp(lnlik+lnprior,eal,ac,score-invp*(b-bp),infaa+invp);
endp;
