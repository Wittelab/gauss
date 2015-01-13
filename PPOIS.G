/* Test program:  ppois(0,1); ppois(1,0); ppois(0,0); end; */
/* This proc computes the exact Poisson probabilities using logfactorials
   Inputs: k =  nonnegative integer
           mu = Poisson mean
   Output: Pr{Poisson(mu) .eq k} */
proc ppois(k,mu);
local test0;
test0 = mu .eq 0;
retp(exp(-mu + k.*ln(mu+test0) - (1-test0).*lnfact(k)) - (test0 .and k .gt 0));
endp;
