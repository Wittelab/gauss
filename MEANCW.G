/* This proc returns the w-weighted column means of the input matrix x.
   w may have either 1 column, for equal weights across columns of x,
   or cols(x) columns, for differing weights across columns of x. */
proc meancw(x,w);
retp(sumc(x.*w)./sumc(w));
endp;