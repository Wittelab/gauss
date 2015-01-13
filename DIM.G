/* This proc prints and returns the dimensions of the matrix x */
proc (0) = dim(x);
"rows & columns: "$+ftocv(rows(x),1,0)$+", "$+ftocv(cols(x),1,0);
retp;
endp;
