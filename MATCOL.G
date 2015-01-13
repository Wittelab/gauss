/* Test program:  let dat[6,3] = 0 0 3  4 1 1  0 0 3  0 0 3  4 1 1  1 0 3;
let rep = 1 1 2 1 8 1; gdat = matcol(dat~rep,1|2|3);
format /rds 5,0; dat~rep; gdat;
output file = temp.out reset; load dat[1633,10] = \gauss\dat\gedi.asc;
t0 = date; gdat = matcol(dat,0);
"Total run time matcol: ";; etstr(ethsec(t0,date)); sumc(gdat[.,cols(gdat)]);
format /rds 5,0; gdat;
end; */
/* This proc sorts and collapses together those rows of a matrix
   that have identical values for all the columns in ccol.
   Values in the remaining columns will be added together.
   Inputs: x = matrix to be collapsed
           ccol = columns to collapse over (0 for all columns)
   Output: xcol = x collapsed and sorted over ccol
                  -- if x is collapsed over all columns, xcol will
                     have a repetition count added as its final column
*/
proc matcol(x,ccol);
local i,acol,xcol,j,k;
if ccol[1] eq 0 or rows(ccol) eq cols(x);
   ccol = seqa(1,1,cols(x)); x = x~ones(rows(x),1); acol = cols(x);
   else; i = seqa(1,1,cols(x)); acol = delif(i,sumr(i .eq ccol'));
endif;
x = sortmc(x,ccol);
xcol = x[1,.];
j = 1; k = 1;
do until j ge rows(x); j = j + 1;
   if prodr(x[j,ccol] .eq x[j-1,ccol]);
      xcol[k,acol] = xcol[k,acol] + x[j,acol];
      else; k = k+1; xcol = xcol|x[j,.];
   endif;
endo;
retp(xcol);
endp;
