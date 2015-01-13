/* This proc returns the negative part of the elements of the matrix x. */
proc negpart(x);
retp((x .lt 0).*x);
endp;

