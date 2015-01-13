/* This procedure generates a vector y of Poisson variates with means mu
   -- truncation is done at far right tail
   -- efficient only if all values of mu are not too large (under 10) */
proc rndpois(mu);
local mumax,kmax,k,cp;
if sumc(mu .lt 0); "ERROR: negative mu passed to rndpois"; end; endif;
/* maximum mu: */ mumax = maxc(mu);
/* maximum value: */ kmax = max(5,ceil(mumax + 4*sqrt(mumax)));
/* truncated range of values: */ k = seqa(0,1,kmax)';
/* expand mu and k to matrices: */
   mu = mu[.,ones(kmax,1)];
   k = k[ones(rows(mu),1),.];
/* cumulative probabilities: */ cp = cdfchic(2*mu,2*(k + 1));
retp(sumr(rndu(rows(mu),1) .ge cp));
endp;
