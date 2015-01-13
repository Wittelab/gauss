/* Test program:
x = seqa(1,1,12); let names = a;
 x = reshape(x,2,6); let names = abcdef bcdef cdef def ef f;
 xint = inter(x); namint = ninter(names);
format /rdn 10,0; x; xint; $namint';
end; */
/* This proc creates interaction (product) terms for regression programs.
   NINTER.G generates names for xint.
   Inputs: x = matrix of terms requiring interactions (products)
   Outputs: xint = matrix of products of x columns (excluding self-products)
   NOTE: If input matrix x has only one column, the proc will return {}.
*/
proc inter(x);
local xint,c,k;
xint = {};  c = cols(x); 
if c eq 1; retp(xint); endif;
k = 1;
do until k eq c; k = k + 1;
   xint = xint~(x[.,(k-1)]*~x[.,k:c]);
endo;
retp(xint);
endp;

