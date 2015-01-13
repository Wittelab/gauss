/* This proc finds MLE of odds ratio correcting for misclassification with known
   classification matrix.
   x = misclassified regressors
   ax = cell counts at x levels
   m = binomial totals (used for binomial only)
   q = classification matrix
   bold = initial b
   bin = 1 for binomial counts, 0 for Poisson
*/
proc (4) = mlmclass(x,ax,m,q,offset,bold,bin);
local nr,odd,iter,b,pt,vt,pm,px,mu,vx,d,bcov,u;
nr = rows(x); odd = seqa(1,2,nr/2);
iter = 0;
if bin; /* binomial */
do until iter ge 100; iter=iter+1;
   pt = expit(offset+x*bold); vt = m.*pt.*(1-pt); pm = q*pt; 
   px = pm; mu = m.*px; vx = mu.*(1-px);
   d = q*(vt.*x); d = d[odd,.]; 
   bcov = invpd(d'((1/vx[odd]).*d)); 
   b = bold + bcov*(d'((ax[odd]-mu[odd])./vx[odd]));
   if prodc(abs(b-bold) .lt .000001); break; endif; bold = b;
endo;
retp(b,bcov,-2*(ax'ln(pm)),iter);
else; /* Poisson */
do until iter ge 100; iter=iter+1;
   vt = exp(offset+x*bold); mu = q*vt;
   d = q*(vt.*x); 
   bcov = invpd(d'((1/mu).*d)); 
   b = bold + bcov*(d'((ax-mu)./mu));
   if prodc(abs(b-bold) .lt .000001); break; endif; bold = b;
endo;
if 0; /* alternate form of OR standard error */
      let u = 1 -1 -1 1; u = u./vt; q = inv(q);
      if sumc(abs(offset)) eq 0; sqrt(u'q*(mu.*q')*u); endif; 
endif;
retp(b,bcov,-2*(ax'ln(mu)-sumc(mu)),iter);
endif;
endp;
end;
