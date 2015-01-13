/* Score and deviance tests bt against maximum for coarsened data
   lnpostp = maximum value of posterior/likelihood
   xr = design columns not tested
   xt = design columns whose coefficients are being tested
   bt = fixed values being tested
   offset (0 if none)
   aa = observed coarsened counts
   cn = numerator coarsening matrix 
        -- maps latent counts to observed counts
   cg = multinomial grouping matrix for coarsened (observed) counts
        -- maps observed counts to fixed margins
   bp, invp = normal-prior mean and covariance (0 for no prior) 
   em = 1 for EM, 0 for Gauss-Newton
   irand = random variables in xr (0 for Poisson)
   bin = initial values for xr coefficients
*/
proc pdevc(lnpostp,xr,xt,bt,offset,aa,cn,cg,bp,invp,em,irand,bin);
local offsets,br,bcovr,ear,acr,mu,lnpostr,score,inf,x2,p; 
offsets = xt*bt+offset;
if em;
    {br,bcovr,ear,acr,lnpostr} = expemc(xr,aa,cn,cg,offsets,irand,bp,invp,bin);
else;
    {br,bcovr,ear,acr,lnpostr} = expgnc(xr,aa,cn,cg,offsets,irand,bp,invp,bin);
endif;
if 0;
   x2 = score'bcovr*score; "score P: ";; cdfchic(x2,rows(bt));;",";;
endif;
x2 = -2*(lnpostr-lnpostp); " dev P  :";; p = cdfchic(x2,rows(bt));
retp(p);
endp;
