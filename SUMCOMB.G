/* Test program:
let x[5,2] = 1 2 3 4 5 6 7 8 9 10 ; x;
print; xc = sumcomb(x,2); xc;
end;
*/
/* This proc recursively constructs the matrix of all possible sums of m1 rows
   taken from the matrix x.
   WARNING: rows(x) must be ge m1, and m1 must be ge 1.
*/
proc sumcomb(x,m1);
/* The trivial cases in which m1 eq 1 (in which case sumcomb = x) or
   m1 eq n (in which case sumcomb = sumc(x)') are disposed of first: */
if m1 eq 1; retp(x); endif;
local n; n = rows(x);
if n eq m1; retp(sumc(x)'); endif;
/* For all other cases, sumcomb is the concatenation of two matrices:
   1) the matrix of all possible sums of m1 rows taken from the first n-1
      rows of x; and
   2) the sum of the last row of x and the matrix of all possible sums of
      m1-1 rows taken from x. */
retp(sumcomb(x[1:(n-1),.],m1)|(sumcomb(x[1:(n-1),.],m1-1) + x[n,.]));
endp;
