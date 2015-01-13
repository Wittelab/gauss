/* This proc creates a diagonal matrix with diagonal vector x. */
proc diagm(x);
retp(diagrv(eye(maxc(rows(x)|cols(x))),x));
endp;