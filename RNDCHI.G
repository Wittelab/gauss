/* This proc generates chi-squared variates by adding gamma and
   Gaussian variates.
   Inputs: df = degrees of freedom
           m,n = dimensions of output matrix
   Output: x2 = m by n matrix of independent chi-squared variates.
*/
proc rndchi(df,m,n);
if df eq 1; retp(rndn(m,n)^2); endif;
if df % 2; retp(reshape(2*sumc(-ln(rndu((df-1)/2,m*n))),m,n) + rndn(m,n)^2);
     else; retp(reshape(2*sumc(-ln(rndu(df/2,m*n))),m,n));
endif;
/* More simply, one could just sum df sqared Gaussian variables:
   retp(reshape(sumc(rndn(df,m*n)^2),m,n));
   but this requires up to twice the number of generated variates. */
endp;
