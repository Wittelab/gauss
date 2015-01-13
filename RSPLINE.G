/* This proc creates regression-spline variables (excluding constant).
   Inputs: x = column vector to convert to spline variables
           k = column vector of knots
           m = order (maximum power) of spline
           rest = 0 for no tail restriction (default), 1 to restrict lower
                    tail to order m-1, 2 to restrict upper tail to order m-1,
                    3 for both tails restricted
           WARNING: USE OF KNOTS OUT OF INPUT-VECTOR RANGE WILL PRODUCE
                    RANK-DEFICIENT OUTPUT MATRICES --
               if minc(x) gt minc(k), collinear columns may be returned;
               if maxc(x) le maxc(k), zero colums will be returned.
           Notes: 1) If m or rest not positive integers,
                     the input matrix will be returned.
                  2) rest will be ignored if m < 2.
   Outputs: matrix of spline variables for x (first column is x).
*/
proc rspline(x,k,m,rest);
if m lt 0 or m ne trunc(m) or rest lt 0 or rest ne trunc(rest);
   "SPLINE.G WARNING: m AND rest MUST BE POSITIVE INTEGERS"; retp(x); endif;
if m eq 0; /* return incremental indicators: */ retp(x .gt k');
   elseif m eq 1; x = x~pospart(x-k');
   elseif m gt 1;
       x = x~(x^seqa(2,1,m-1)')~(pospart(x-k')^m);
       if rest ge 2;
          x = x[.,1:(m-1)]~(x[.,m:(cols(x)-1)]-x[.,cols(x)]); endif;
       if rest eq 1 or rest ge 3; x = x[.,1:(m-1) (m+1):cols(x)]; endif;
   else; "ERROR IN SPLINE.G: Invalid order for spline"; retp(x);
endif;
retp(x);
endp;
