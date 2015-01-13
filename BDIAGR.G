/* This proc creates a block diagonal matrix from a matrix
   that is a row of square blocks.
   Input: x = a matrix of the form x1~...~xn, where the xi
              are k-by-k square matrices
   Output: d = an nk-by-nk square block diagonal matrix with i'th diagonal
               block xi.
*/
proc bdiagr(x);
local n,k,d,i,s1,sk;
k = rows(x);
n = cols(x)/k;
if n*k ne cols(x);
   "ERROR IN BDIAGR.G: cols(x) is not divisible by rows(x)"; end;
endif;
d = zeros(n*k,n*k); i = 0;
do until i ge n; i = i + 1;
   s1 = (i-1)*k + 1; sk = i*k;
   d[s1:sk,s1:sk] = x[.,s1:sk];
endo;
retp(d);
endp;
