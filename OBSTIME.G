/* This proc returns the time observed in each interval on each subject, where
   Inputs: te = starts of observation
           tx = ends of observation
           lcutp = start of first interval or starts of all intervals
           ucutp = ends of all intervals or end of last interval
   Output: pt = person-time observed in each interval
   NOTES: 1) pt has length = no. of subjects times number of intervals.
          2) THERE MUST BE AT LEAST TWO INTERVALS
*/
proc obstime(te,tx,lcutp,ucutp);
if rows(te) eq 1; te = te*rows(tx); endif;
if rows(lcutp) eq 1; lcutp = lcutp|ucutp[1:(rows(ucutp)-1)]; endif;
if rows(ucutp) eq 1; ucutp = lcutp[2:rows(lcutp)]; endif;
retp(pospart(minr((rows(tx).*.ucutp)~(tx.*.rows(ucutp))
 - maxr((rows(te).*.lcutp)~(te.*.rows(lcutp))));
endp;
