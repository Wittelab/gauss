/* This procedure finds the medians of the columns of x.
   Runs in about 76% the time of median(x). */
proc medianc(x);
local im,evenr,i,m,v;
if rows(x) eq 1; retp(x'); endif;
im = rows(x)/2;
/* indicate even no. of rows: */ evenr = im eq floor(im);
/* im is middle row index for odd no. rows, lower of middle two indices
   for even no. rows: */ im = im + (1 - evenr)/2;
i = 0; m = zeros(cols(x),1);
do until i ge cols(x); i = i + 1;
   v = sortc(x[.,i],1);
   m[i] = ((v[im] + v[im+evenr])/2);
endo;
retp(m);
endp;
