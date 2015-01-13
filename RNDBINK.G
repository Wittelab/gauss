/* This procedure generates a k-vector y of binomial
   variates with probability p and total n */
proc rndbink(p,n,k);
retp(sumc(rndu(n,k) .lt p));
endp;
