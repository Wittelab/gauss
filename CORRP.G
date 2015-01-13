/* This proc returns the matrix of all pairwise correlations between columns
   of x and columns of y.
   output:  correlation matrix with columns of x corresponding to rows
            and columns of y corresponding to columns; thus, the returned
            matrix will have cols(x) rows and cols(y) columns.
*/
/* Test program:
let x[3,2] = 1 1  2 3  3 2; let y[3,3] = 1 1 3  3 2 2  2 3 1;
format 10,3; corrp(x,y); end; */
proc corrp(x,y);
local r; r = corrx(x~y);
retp(r[1:cols(x),cols(x)+1:cols(r)]);
endp;
