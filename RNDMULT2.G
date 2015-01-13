/* This program generates a multinomial variate
   with distribution p and total n.
   Faster than rndmultn.g when n > 55. */
proc rndmult2(p,n);
retp(counteq((sumc(rndu(1,n) .le cumsumc(rev(p)))),seqa(1,1,rows(p)),0));
endp;

