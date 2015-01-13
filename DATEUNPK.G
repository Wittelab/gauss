/* Test program: let d = 030709 150987 112089; dateunpk(d); end; */
/* This proc unpacks a column vector d of dates to a 3-col matrix.
   Input: d = dates
   Output: date matrix with mo-day-yr columns if d has mmddyy format,
           day-mo-yr columns if d has ddmmyy format.
*/
proc dateunpk(d);
retp(trunc((d/10000)~((d%10000)/100)~(d%100)));
endp;
