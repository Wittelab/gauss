/* For each element of y, this proc counts the number of occurrences in x
   that are > or = the element. Equivalent to sumr(x'.ge y).
   Inputs: x,y, and
           a = repetition count for x (set to 0 if no rep count)
   Output: no. of occurrences in x that are ge y
   NOTE: Output has one row for each row of y.
         y MUST be sorted in ascending order.
 */
proc countge(x,y,a);
if sumc(a);
   retp(rev(cumsumc(countwts(-x,-rev(y),a))));
else;
   retp(rev(cumsumc(counts(-x,-rev(y)))));
endif;
endp;
