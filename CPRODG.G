/* Test program:  let dat[5,3] = 2 1 1  4 2 1  2 2 3  3 5 3  2.5 1.5 3;
gdat = cprodg(dat,3); format /rds 5,1; dat~gdat; end; */
/* This proc creates a matrix of within-group cumulative products.
   Groups are defined by those rows of the input matrix
   that have identical values for all the columns in cg.
   NOTE: x must be sorted on cg overall and as desired within levels of cg
   Inputs: x = matrix to be operated on
           cg = columns that identify the groups for cumulative products
   Output: xp = x with cg-specific cumulative products replacing columns
                with indices not in cg
*/
proc cprodg(x,cg);
local i,wg,xp,j;
i = seqa(1,1,cols(x)); wg = delif(i,sumr(i .eq cg'));
xp = x;
j = 1;
do until j ge rows(x); j = j + 1;
   if prodr(x[j,cg] .eq x[j-1,cg]);
      xp[j,wg] = x[j,wg].*xp[j-1,wg];
      else; xp[j,wg] = x[j,wg];
   endif;
endo;
retp(xp);
endp;
