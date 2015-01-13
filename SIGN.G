/* This proc returns -1 for negative arguments, 1 for positive arguments,
   0 for zero arguments. */
proc sign(x);
retp((x .gt 0) - (x .lt 0));
endp;
