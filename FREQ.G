/* This proc does frequency counts on the columns of x. Column
   proportions will also be printed if only 1 count per column is specified,
   Inputs: x = data matrix.
           rep = column index(es) of repetition count(s)  (0 if none);
                 if multiple rep counts are specified, all will be printed.
           names = names of columns of x (0 if no printed output desired).
   Outputs: nvals = no. of distinct values found in each column of x.
            uvals = vector containing the successive subvectors of unique
                    values for each column of x.
            freqs = vector containing the successive subvectors of
                    frequency counts for each column of x;
                    will be matrix if multiple rep counts are specified.
            is = starts of code & frequency subvectors in uvals & freqs.
            ie = ends of code & frequency subvectors in uvals & freqs.
*/
proc (5) = freq(x,rep,names);
local uvals,nvals,freqs,j,u,f,k,ie,mask,fmt;
uvals = {}; nvals = {}; freqs = {}; j = 0;
do until j ge cols(x); j = j + 1;
   u = unique(x[.,j],1);
   uvals = uvals|u; nvals = nvals|rows(u);
   if rows(rep) eq 1;
      f = counts(x[.,j],u);
      else; k = 0; f = {};
        do until k ge cols(rep); k = k + 1;
           f = f~countwts(x[.,j],u,rep[.,k]);
        endo;
   endif;
   freqs = freqs|f;
   if 0$+names[1]; format /rds 10,0; $names[j];;" has";;nvals[j];;
      let fmt[3,3] = "*.*lf " 10 3 "*.*lf " 10 0 "*.*lf" 9 3;
      if cols(f) eq 1; " values with counts & proportions:";
               call printfm(u~f~(f/sumc(f)),1~1~1,fmt);
         else; " values with counts:";
               mask = ones(1,cols(f)); fmt[2,2] = 5;
               call printfm(u~f,1~mask,fmt[1~(2*mask),.]); format 5,0;
      endif; print;
   endif;
endo;
if 0$+names[1]; "Total:     ";; sumc(f)'; endif;
ie = cumsumc(nvals);
retp(nvals,uvals,freqs,ie-nvals+1,ie);
endp;
