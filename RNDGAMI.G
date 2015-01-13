/* This program generates gamma variates for integral values of the shape
   parameter and unit values of the scale parameter by summing exponential
   variates. Not recommended for large values of a. See Devroye, 1986, p. 405.
   NOTE: To generate a gamma(a + 1/2,1) variate, add a chi-squared(1)/2 variate
         to a gamma(a,1) variate, or use a chi-squared(2a+1)/2 variate.
   Inputs: a = integral value of scale parameter
           m,n = rows & columns for output.
   Output: x = m-by-n matrix of gamma(a,1) variates.
*/
proc rndgami(a,m,n);
retp(reshape(sumc(-ln(rndu(a,m*n))),m,n));
endp;

