/* ln-posterior evaluation with coarsened loglinear Poisson counts.
   Inputs:
   aa = coarsened counts
   x = complete-data design matrix
   cn = coarsening matrix
   b = parameter value at which to evaluate posterior
   bp, invp = normal-prior mean and inverse covariance (0 if none)
   Outputs: 
   ln posterior (ln likelihood if no prior)
   eal = expected complete (latent) data 
   ac = expected complete data given observed aa
   ln posterior derivative (score if no prior)
   expected negative ln posterior Hessian (expected information if no prior)
*/
proc (5) = postpc(aa,x,cn,b,offset,bp,invp);
local eal,eaa,lnlik,d,score,inf,ac,lnprior;
   /* complete-data expectations: */ eal = exp(offset+x*b);
   /* coarsened expectations: */ eaa = cn*eal; 
   /* ln likelihood: */ lnlik = aa'ln(eaa) - sumc(eaa+lnfact(aa));
   /* derivative of ea: */ d = cn*(eal.*x);
   /* lnlik gradient: */ score = d'(aa./eaa - 1); 
   /* expected information: */ inf = d'(d./eaa);
   ac = (cn'(aa./(cn*eal))).*eal; 
   lnprior = -((b-bp)'invp*(b-bp))/2;
   retp(lnlik+lnprior,eal,ac,score-invp*(b-bp),inf+invp);
endp;
