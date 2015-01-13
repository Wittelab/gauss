/* This proc partially sums the columns of matrix x
   according to the sequence given in m.
   Inputs: x = matrix to be partially summed
           m = vector giving the number of rows of x contributing to
               each partial sum, in sequence
   Outputs: z = matrix of partial sums
   NOTE: rows(z) = rows(m) */
proc mpsumc(x,m);
local t,s,i,z;
t = cumsumc(m); s = t - m + 1;
z = {}; i = 0;
do until i ge rows(m);  i = i + 1;
   z = z|sumc(x[s[i]:t[i],.])';
endo;
retp(z);
endp;
