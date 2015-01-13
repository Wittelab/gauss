/* This proc checks that range spanned by each column j of x is within
   the range of the corresponding column of r (or the single column of r).
   Inputs: x = matrix to be checked
           r = range matrix
           names = names of x columns (0 if no printed output desired)
   Outputs: xflag = 2-by-cols(x) matrix, with columns corresponding to
                   x columns as follows:
                   row 1: 1 indicates lower range violation, 0 otherwise
                   row 2: 1 indicates upper range violation, 0 otherwise
            xrange = 2-by-cols(x) matrix with minc(x) in row 1 and maxc(x)
                    in row 2.
*/
proc (2) = rangechk(x,r,names);
local xrange,xflag,nv;
xrange = (minc(x)~maxc(x))';
xflag = (xrange[1,.] .lt minc(r)')|(xrange[2,.] .gt maxc(r)');
if 0$+names[1]; /* report range violations: */ nv = sumc(xflag');
   if nv[1]; print;
      " lower range violation(s) for "; $selif(names,xflag[1,.]')';
      " which have minima "; selif(xrange[1,.]',xflag[1,.]')';
   endif;
   if nv[2];
      " upper range violation(s) for "; $selif(names,xflag[2,.]')';
      " which have maxima "; selif(xrange[2,.]',xflag[2,.]')';
   endif;
endif;
retp(xflag,xrange);
endp;
