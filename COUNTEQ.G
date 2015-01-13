/* For each element of y, this proc counts the number of occurrences in x
   that equal the element. Equivalent to sumr(x'.eq y).
   Inputs: x,y, and
           a = repetition count for x (set to 0 if no rep count)
   Output: no. of occurrences in x that eq y
   NOTES: Output has one row for each row of y.
          y MUST be sorted in ascending order.
 */
proc counteq(x,y,a);
if sumc(a);
   retp(countge(x,y,a) + countle(x,y,a) - sumc(a));
else;
   retp(countge(x,y,0) + countle(x,y,0) - rows(x));
endif;
endp;


