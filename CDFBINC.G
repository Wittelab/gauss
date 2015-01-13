/* This proc computes the exact upper tail of a binomial(p,n) distribution
   using the F distribution (see, e.g., Ling, Amer Statist 1992;46:53-54.
   Inputs: k = upper cut-off point (integer between 0 and n);
           p = Bernoulli probability
           n = no. of Bernoulli trials
   Output: Pr{bin(p,n) .ge k} */
/* Binomial test program: */
/* format 12,8; let n = 3 4 5 5; let k = 2 2 0 5; let p = 0.2;
n'; k'; u = cdfbinc(k,p,n)'; l = cdfbin(k,p,n)'; p = pbinf(k,p,n)';
u; l; u+l-p; p; end; */
proc cdfbinc(k,p,n);
local test1;
test1 = (k .eq 0) .or (p .eq 1);
if sumall(k .gt n) or sumall(k .lt 0) or sumall(p .gt 1) or sumall(p .lt 0);
   "ERROR: p or k out of range in cdfbinc"; end;
   else;
   retp(1 - (1-test1).*cdffc((n-k+1).*p./(k.*(1-p) + test1),2*k,2*(n-k+1)));
endif;
endp;
