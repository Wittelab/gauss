/* This procedure computes the trace of a matrix */
proc tracemat(x);
retp(sumc(diag(x)));
endp;

