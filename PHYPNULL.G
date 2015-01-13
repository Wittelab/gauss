/* This proc computes the null hypergeometric probabilities of
   the 2x2 tables that are the rows of t: */
proc phypnull(t);
local m;
/* margins: */ m = (t[.,1]+t[.,2 3])~(t[.,2 3]+t[.,4]);
retp(exp(sumr(lnfact(m)) - lnfact(sumr(t)) - sumr(lnfact(t))));
endp;
