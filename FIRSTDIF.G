/* This proc returns the matrix of first differences for the columns
   of x. Note that the returned matrix has one less row than x. */
proc firstdif(x);
retp(x[1:rows(x)-1,.]-x[2:rows(x),.]);
endp;

