/* This proc computes the exact upper tail of a Poisson distribution
   using the chi-squared distribution (see, e.g., Brownlee)
   Inputs: k = upper cut-off point (nonnegative integer);
           mu = Poisson mean
   Output: Pr{Poisson(mu) .ge k} */
proc cdfpoisc(k,mu);
local test0;
test0 = (k .eq 0);
if sumall(k .lt 0) or sumall(mu .lt 0);
   "ERROR: Negative k or mu passed to cdfpoisc"; end;
   else; retp(1 - (1-test0).*cdfchic(2*mu,2*(k + test0)));
endif;
endp;

