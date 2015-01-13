/* Test program: x = 100*rndu(11,2); c = pctile(x,3|4);
 format /rds 9,3; c; sortc(x[.,1],1)~sortc(x[.,2],1); end; */
/* This proc returns the boundaries for percentile (quantile) categories,
   including the maximum as the top boundary.
   WARNING: There is no provision for ties.
   NOTES: 1) If rows(x) is not divisible by ncat, the remaining units are
             allocated one-by-one to categories starting from the highest.
          2) Boundaries are included in the category below them.
   Inputs: x = matrix whose columns are to be categorized
           ncats = no. of categories for each x column (e.g., 4 for quartiles)
                  -- if vector, rows(ncats) must equal cols(x)
   Output: bounds = category boundaries: bounds[1:ncats[i],i]
                  (which is the first ncats[i] rows of column i of bounds)
                  contains the boundaries for column i of x; any
                  remaining rows of column i of bounds are padded with zeros.
*/
proc pctile(x,ncats);
local n,bounds,i,nc,r,s;
n = rows(x); if rows(ncats) eq 1; ncats = ncats*ones(cols(x),1); endif;
bounds = zeros(maxc(maxc(ncats)),cols(x)); i = 0;
do until i ge cols(x); i = i + 1; nc = ncats[i];
   /* remainder: */ r = n%nc;
   if r; s = zeros(nc-r,1)|ones(r,1); else; s = 0; endif;
   s = cumsumc(floor(n/nc)*ones(nc,1) + s);
   x[.,i] = sortc(x[.,i],1);
   bounds[1:nc,i] = x[s,i];
endo;
retp(bounds);
endp;
