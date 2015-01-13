/* Test program:  format /rdn 6,0; vecomb(seqa(1,1,5),3); end; */
/* This proc recursively constructs the matrix of all possible row vectors of
   length m1 with elements taken from the column vector x, in matrix order
   (run test program to see).
   NOTE: To save computation, vecomb.g does NOT check that n ge m1 ge 1,
         required.
*/
proc vecomb(x,m1);
if m1 eq 1; retp(x); endif;
local n; n = rows(x);
if n eq m1; retp(x'); endif;
local c; c = combin(n-1,m1-1);
retp(((x[1]*ones(c,1))~vecomb(x[2:n],m1-1))|vecomb(x[2:n],m1));
endp;
