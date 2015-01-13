/* This proc prints frequency counts on the columns of x. 
   Inputs: x = data matrix.
           rep = column index(es) of repetition count(s)  (0 if none);
                 if multiple rep counts are specified, all will be printed.
           names = names of columns of x (0 if no printed output desired).
   Outputs: nvals = vector of number of distinct values found in each column of x.         
   FORMAT CONTROLLED BY CALLING PROGRAM!
*/
proc freqmin(x,rep,names);
local nvals,i,u;
if rows(rep) eq 1; rep = ones(rows(x),1); 
   elseif rows(rep) ne rows(x); "FREQMIN: ERROR IN rep"; end; endif;
"VALUES AND FREQUENCIES"; nvals = zeros(cols(x),1); i=0;
do until i ge cols(x); i = i + 1;
   u=unique(x[.,i],1)'; nvals[i] = cols(u);
   $(""$+names[i]$+" has "$+ftocv(nvals[i],1,0)$+" values --"); 
   " values: ";; u; " counts: ";; rep'(x[.,i] .eq u);
endo;
retp(nvals);
endp;
