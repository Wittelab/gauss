/* For each element of y, this proc counts the number of occurrences in x
   that are < or = the element. Equivalent to sumr(x'.le y).
   Inputs: x,y, and
           a = repetition count for x (set to 0 if no rep count)
   Output: no. of occurrences in x that are le y
   NOTE: Output has one row for each row of y.
         y MUST be sorted in ascending order.
 */
proc countle(x,y,a);
if sumc(a);
   retp(cumsumc(countwts(x,y,a)));
else;
   retp(cumsumc(counts(x,y)));
endif;
endp;

