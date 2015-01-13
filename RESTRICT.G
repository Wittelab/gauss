/* This proc restricts (truncates) a vector x by lb below and ub above.
   lb and ub may be scalar or vector of same dimension as x. */
proc restrict(x,lb,ub);
retp(min(max(x,lb),ub));
endp;
