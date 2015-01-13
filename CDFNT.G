/* This proc computes the two-tailed area for the standard normal distribution.
   Input: z = standard normal variate.
   Output: area in tails of standard normal distribution outside of +/- z */
proc cdfnt(z);
retp(2*cdfnc(abs(z)));
endp;