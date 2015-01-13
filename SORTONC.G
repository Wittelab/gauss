/* Test program:  n = 20; x = ceil(10*rndu(2*n,1));
v = permuter(unique(x,1)); v = v[1:5]|100|v[6:10];
xs = sortonc(x,1,v); format /rds 5,0;
x[1:n]~xs[1:n]~x[n+1:2*n]~xs[n+1:2*n]; v'; end; */
/* This proc sorts the rows of the numeric matrix x based on the ordering of
  values for x[.,j] given in the vector v; v MUST include all the distinct
  values in x[.,j] with no duplications or else the program will terminate
  with an error message.
  Inputs: x = matrix to be sorted
          j = column of x to sort on
          v = list of all distinct values in x, in the order desired for
              the sorted matrix (top to bottom rows)
  Output: x = sorted version of input x
*/
proc sortonc(x,j,v);
local ux,uv;
ux = unique(x[.,j],1); uv = unique(v,1);
if rows(v) ne rows(uv);
   "ERROR IN SORTONC.G: v has duplicate values."; end; endif;
if sumall(ux .eq uv') ne rows(ux);
   "ERROR IN SORTONC.G: v does not contain all the values of x[.,j]."; end;
endif;
retp(x[sortind(sumr((x[.,j] .eq v').*seqa(1,1,rows(v))')),.]);
endp;
