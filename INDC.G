/* This proc takes a matrix x and returns a matrix
   xind with the same number of rows, each column an indicator for whether
   the x value in that row of col k equals the corresponding value in codes.
   NINDC.G generates names for the columns of xind.
   Inputs: x = categorical data matrix
           codes = vector of possible values for x components
           s = vector of indices at which the set of codes for column
               i of x begin:
           Ordinarily, rows(s) = cols(x). If s = 1 and cols(x) > 0,
           the same values (those in codes) will be used for each x column
   Output: xind = indicator matrix
   Note:  cols(xind) = rows(codes)
*/
proc indc(x,codes,s);
local k,e,xind;
if rows(s) eq 1 and cols(x) gt 1; /* expand s and codes: */
   s = seqa(1,rows(codes),cols(x)); codes = ones(cols(x),1).*.codes;
endif;
if rows(s) gt 1; e = (s[2:rows(s)] - 1)|rows(codes);
   else; e = rows(codes); endif;
k = 0; xind = {};
do until k ge cols(x); k = k + 1;
   xind = xind~(x[.,k] .eq codes[s[k]:e[k]]');
endo;
retp(xind);
endp;
