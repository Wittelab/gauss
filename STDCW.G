/* This proc returns the w-weighted column standard deviations
   of the input matrix x.
   w may have either 1 column, for equal weights across columns of x,
   or cols(x) columns, for differing weights across columns of x. */
proc stdcw(x,w);
retp(sqrt((w'(x-meancw(x,w)')^2)./sumc(w)));
endp;