/* This proc takes a matrix x and returns an incremental-indicator matrix
   xind with the same number of rows as x, and with a block of columns
   corresponding to each column of x. Each column in block k of xind indicates
   whether the x[.,k] value in that row is greater than or equal to a unique
   observed value of x[.,k]. The indicators for a column are ordered to
   correspond to ascending values in the column.
   NINDC2.G will generate names for the columns of xind.
   Inputs: x = categorical data matrix
   Output: xind = incremental indicator matrix
*/
/* Test program:
let x[5,2] = 1 0  2 1  2 2  3 2  3 1;
xi = indinc(x); x[.,1]~xi[.,1 2]; x[.,2]~xi[.,3 4]; end; */
proc indinc(x);
local k,u,xind;
k = 0; xind = {};
do until k ge cols(x); k = k + 1;
   u = unique(x[.,k],1);
   if rows(u) eq 1; "COLUMN ";; k;; " HAS ONLY ONE UNIQUE VALUE."; end; endif;
   xind = xind~(x[.,k] .ge u[2:rows(u)]');
endo;
retp(xind);
endp;


