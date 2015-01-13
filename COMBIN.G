/* Computes the combinatorial coefficients n!./(k!.*(n-k)!) */
proc combin(n,k);
retp(round(exp(lnfact(n) - lnfact(k) - lnfact(n-k))));
endp;

