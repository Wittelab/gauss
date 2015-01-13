/* Test program:
let s[3,2] = 1 2  2 3  2 1;
x = reshape(seqa(0,1,9),3,3); format /rdn 10,0; x; s;
seleltc(reshape(seqa(0,1,9),3,3),s);
end; */
/* This proc selects the column elements indicated in s from x.
   s must have the same number of rows as x,
   and the values in s must be from 1,2,...,cols(x).
   Output: sx = matrix of same dimension of s with sx[i,j] = x[i,s[i,j]].
*/
proc seleltc(x,s);
local ic,j,sx;
if (minc(minc(s)) lt 1) or (maxc(maxc(s)) gt cols(x));
   "ERROR IN SELELTC.G: element of s out of range -- "; format /rdn 10,3;
   " Min & max elements of s are ";; minc(minc(s));; maxc(maxc(s)); end;
endif;
ic = seqa(1,1,cols(x))';
j = 0; sx = {};
do until j ge cols(s); j = j + 1;
   sx = sx~sumr(x.*(s[.,j] .eq ic));
endo;
retp(sx);
endp;

