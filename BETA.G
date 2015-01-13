/* computes beta function */
proc beta(a,b);
retp(exp(lngamma(a)+lngamma(b)-lngamma(a+b)));
endp;