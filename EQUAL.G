/* This proc checks equality of x and y */
proc equal(x,y);
local n; n = 0; 
"X: ";; dim(x); "Y: ";; dim(y);
if rows(x) ne rows(y) or cols(x) ne cols(y); "dim(X) not equal dim(Y)"; retp(-1); endif;
n = sumall(x .ne y); "No. X not = Y: "$+ftocv(n,1,0);
retp(n);
endp;