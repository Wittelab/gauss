/* Test program:
x = seqa(1,1,4); m = rev(x); x = x~(x+4); mexpand(x,m); end; */
/* This proc expands the matrix x by creating a matrix with row k of x
   repeated rep[k] times
   Inputs: x = matrix to be expanded
           rep = vector of repetition counts for each row of x
           NOTE: All elements of m must be positive.
   Output: expanded matrix
*/
proc mexpand(x,rep);
local ix,i,s;
ix = zeros(sumc(rep),1); i = 0; s = 0;
do until i ge rows(x); i = i + 1;
   ix[s+1:s+rep[i]] = (i*ones(rep[i],1)); s = s + rep[i];
endo;
retp(x[ix,.]);
endp;

