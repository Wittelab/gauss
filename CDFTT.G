/* This proc computes the two-tailed area for a t distribution
   Input: t = t variates
   Output: area in tails of t distribution outside of +/- t */
proc cdftt(t,df);
retp(2*cdftc(abs(t),df));
endp;