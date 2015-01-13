/* This proc computes the exact lower tail of a binomial(p,n) distribution
   using the F distribution (see, e.g., Ling, Amer Statist 1992;46:53-54.
   Inputs: k = upper cut-off point (integer between 0 and n);
           p = Bernoulli probability
           n = no. of Bernoulli trials
   Output: Pr{bin(p,n) .le k} */
proc cdfbin(k,p,n);
local test1;
test1 = (k .eq n) .or (p .eq 0);
if sumall(k .gt n) or sumall(k .lt 0) or sumall(p .gt 1) or sumall(p .lt 0);
   "ERROR: p or k out of range in cdfbin"; end;
   else;
   retp(1-(1-test1).*cdffc((k+1).*(1-p)./((n-k).*p + test1),2*(n-k),2*(k+1)));
endif;
endp;
