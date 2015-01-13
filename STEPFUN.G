/* Test program: b = seqa(0,1,10); h = 10*rndu(9,2);
   { x,y } = stepfun(b,h); library pgraph; xy(x,y); end; */
/* This proc returns a graphable step function from the inputs.
   Inputs: b = category bounndaries (including min and max)
           h = heights of steps in categories
   Ouputs: x = sequence(s) of abscissa values
           y = sequences of ordinate values
   NOTE: 1) b must have one more row than h.
         2) x and y will have 2*rows(b) rows.
         3) b and/or y may be matrices, in which case x and/or y
            will be matrices. If both b and y are matrices, they must have
            the same number of columns.
*/
proc (2) = stepfun(b,h);
local s; s = seqa(2,1,rows(b)-2).*.(1|1);
retp(b[1|s|rows(b),.],h[1|1|s,.]);
endp;

