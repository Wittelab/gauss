/* This proc counts the asymmetries in a square matrix x.
   Output = no. of asymmetries, -1 if x not square */
proc asymetry(x);
if rows(x) ne cols(x);
   "ASYMETRY.G: Input matrix not square."; retp(-1); endif;
retp(sumall(x .ne x'));
endp;