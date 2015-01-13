/* This proc takes a matrix x and returns a matrix xc of the same size whose
   elements are the number of cutpoints in cutp that the corresponding element
   of x is greater than.
   Inputs: x = matrix of variables to be categorized
           cutp = vector of cutpoints
           s:  if s is a scalar, it will be assumed that cutp contains only one
               set of cutpoints, to be used for all columns of x;
               otherwise, s should be the vector of indices at which the set of
               cutpoints for column i of x begin, in which case rows(s) should
               equal cols(x).
   Output: xc = matrix of categorized values of x.
   Notes: 1) The possible values for xc[.,i] range from zero to the number of
             cutpoints for x[.,i].
          2) unique(xc[,.i],1) will provide the distinct values of xc[.,i].
*/
proc categ(x,cutp,s);
local i,e,xc;
if rows(s) gt 1;
   e = (s[2:rows(s)] - 1)|rows(cutp);
   i = 0; xc = {};
   do until i ge cols(x); i = i + 1;
      xc = xc~sumr(x[.,i] .gt cutp[s[i]:e[i]]');
   endo;
   else; xc = reshape(sumr(vec(x) .gt cutp'),rows(x),cols(x));
endif;
retp(xc);
endp;
