/* This proc takes a matrix x and returns a matrix
   xind with the same number of rows, each column an indicator for whether
   the x value in that row of col k equals a particular value seen in that
   column. The indicators for a column are ordered to correspond to
   ascending values in the column.
   NINDC2.G generates names for the columns of xind.
   Inputs: x = categorical data matrix
           ref = vector of reference values for x columns -
   NOTE:   No indicator will be generated for the values in ref.
           If ref is scalar, the same ref value will be used for all columns.
   Output: xind = indicator matrix
*/
/* Test program:
let x[5,2] = 1 0  2 1  2 2  3 2  3 1; let ref = 1 0;
xi = indc2(x,ref); x[.,1]~xi[.,1 2]; x[.,2]~xi[.,3 4]; end; */
proc indc2(x,ref);
local k,u,xind;
k = 0; xind = {};
if rows(ref) eq 1; ref = ref*ones(cols(x),1); endif;
do until k ge cols(x); k = k + 1;
   u = unique(x[.,k],1); u = delif(u,u .eq ref[k]);
   xind = xind~(x[.,k] .eq u');
endo;
retp(xind);
endp;
