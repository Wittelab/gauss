/* This proc computes the exact lower tail of a Poisson distribution
   using the chi-squared distribution (see, e.g., Ulm, AJE 1990;131:373-5)
   Inputs: k = lower cut-off point (nonnegative integer);
           mu = Poisson mean
   Output: Pr{Poisson(mu) .le k} */
proc cdfpois(k,mu);
if sumall(k .lt 0) or sumall(mu .lt 0);
   "ERROR: Negative k or mu passed to cdfpoisc"; end;
   else; retp(cdfchic(2*mu,2*(k + 1)));
endif;
endp;


