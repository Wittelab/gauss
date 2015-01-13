/* This proc returns a column vector r
   containing the ranges of each column of the input matrix x */
proc rangec(x);
retp(maxc(x) - minc(x));
endp;