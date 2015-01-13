/* This proc does frequency counts & group means on the columns of x.
   Inputs: x = data matrix
           y = variable for within-group means
           rep = repetition count (0 if none)
           xnames = names of columns of x (0 if no printed output desired)
           yname = name of y (may be 0 if no printed output desired)
   Outputs: nvals = no. of distinct values found in each column of x
            uvals = vector containing the successive subvectors of unique
                    values for each column of x
            freqs = vector containing the successive subvectors of
                    frequency counts for each column of x
            means = vector containing the successive subvectors of
                    subgroup means of y for each column of x
            sds = vector containing the successive subvectors of
                  subgroup std dev of y for each column of x
            is = starts of code & frequency subvectors in uvals & freqs
            ie = ends of code & frequency subvectors in uvals & freqs
*/
proc (7) = freqmean(x,y,rep,xnames,yname);
/* precision of x values: */ declare matrix _xprec = 3;
/* precision of y stats: */ declare matrix _yprec = 3;
local y2,fmt,tot,uvals,nvals,freqs,means,sds,j,u,f,m,s,ie;
y2 = y.*y;
if 0$+xnames[1];
   "VALUES, COUNTS, & PROPORTIONS OF TOTAL,"; format 10,3;
   "AND SUBGROUP MEANS & SAMPLE STANDARD DEVIATIONS FOR "$+yname$+":";
   let fmt[5,3] = "*.*lf" 10 3 "*.*lf" 10 0
                  "*.*lf" 10 3 "*.*lf" 10 3 "*.*lf" 10 3;
   fmt[1,3] = _xprec; fmt[4,3] = _yprec; fmt[5,3] = _yprec;
endif;
if rows(rep) eq rows(x);
   y = rep.*y; y2 = rep.*y2; tot = sumc(rep); else; tot = rows(x); endif;
uvals = {}; nvals = {}; freqs = {}; means = {}; sds = {}; j = 0;
do until j ge cols(x); j = j + 1;
   u = unique(x[.,j],1);
   uvals = uvals|u; nvals = nvals|rows(u);
   if rows(rep) eq 1;
      f = counts(x[.,j],u);
      else; f = countwts(x[.,j],u,rep);
   endif;
   m = countwts(x[.,j],u,y)./f;
   s = (countwts(x[.,j],u,y2) - f.*m.*m)./(f - (f .gt 1));
   freqs = freqs|f;
   means = means|m;
   sds = sds|s;
   if 0$+xnames[1];
      format 10,0; $xnames[j];;" has";;nvals[j];;" distinct values:";
      printfm(u~f~(f/tot)~m~s,ones(1,5),fmt);
   endif;
endo;
if 0$+xnames[1];
   "Total count & mean:";; printfm(tot~(sumc(y)/tot),1~1,fmt[2 3,.]); endif;
ie = cumsumc(nvals);
retp(nvals,uvals,freqs,means,sds,ie-nvals+1,ie);
endp;

