/* This proc computes the exact probability of a binomial(p,n) outcome k
   using the F distribution (see, e.g., Ling, Amer Statist 1992;46:53-54)
   Inputs: k = Binomial outcome (integer between 0 and n);
           p = Bernoulli probability
           n = no. of Bernoulli trials
   Output: Pr{bin(p,n) .eq k} */
proc pbinf(k,p,n);
retp(cdfbin(k,p,n) + cdfbin(n-k,1-p,n) - 1);
endp;
