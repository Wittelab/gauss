/* This procedure generates a matrix of Bernoulli
   variates with probabilities p */
proc rndbern(p);
retp(rndu(rows(p),cols(p)) .le p);
endp;
