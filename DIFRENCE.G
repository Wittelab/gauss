/* This proc computes vertical (row) differences of any order m.
   Inputs: x = matrix to be differenced
           m = order of differencing
    Output: x = difference matrix (differences of x rows)
                (this will have m fewer rows than the input x) */
proc difrence(x,m);
local i; i = 0;
do until i ge m; i = i+1;
   x = x[1:rows(x)-1,.] - x[2:rows(x)];
endo;
retp(x);
endp;

