/* This proc computes the cumulative products of rows of a matrix. */
proc cumprodr(x);
retp(cumprodc(x')');
endp;
