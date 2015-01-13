/* This proc computes the cumulative sums of rows of a matrix. */
proc cumsumr(x);
retp(cumsumc(x')');
endp;
