/* This proc creates a block diagonal matrix from a matrix that is a column of
  blocks, which need not be square. See bdiagc.g or bdiagr.g for square blocks.
   Input: x = a matrix of the form x1|...|xn, where the xi
              are matrices of constant dimension
          k = rows(xi)
   Output: d = block diagonal matrix with i'th diagonal block xi
               rows(d) = rows(x), cols(d) = (r/k)*cols(x)
*/
proc bdiag(x,k);
local r,c,n,d,i,r1,rk,c1,cc;
r = rows(x); c = cols(x);
n = r/k;
if n*k ne r;
   "ERROR IN BDIAGM.G: rows(x) is not divisible by k"; end;
endif;
d = zeros(r,c*n); i = 0;
do until i ge n; i = i + 1;
   r1 = (i-1)*k + 1; rk = i*k; c1 = (i-1)*c + 1; cc = i*c;
   d[r1:rk,c1:cc] = x[r1:rk,.];
endo;
retp(d);
endp;
