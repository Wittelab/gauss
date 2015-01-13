/* This program creates a crosstabulation matrix from a data matrix, and
   prints the matrix if variable names are supplied.
   Inputs: x = matrix of variables to be crosstabulated
           ncat = vector giving number of code values to be used for each
                  variable -- if 0, the program will use all unique values
                  in each x column as the codes.
           codes = vector of length sumc(ncat) giving in succession the
                   code values for each variable (will be ignored if ncat = 0)
           rep = repetition count (0 if none) --
                 may have multiple columns, in which case t will have
                 multiple count columns following the code columns.
           xnames = names of x columns (0 if no printed output desired)
           WARNING: PRINT FORMAT IS CONTROLLED BY CALLING PROGRAM!
   Outputs: t = prodc(ncat)-row by (cols(x)+cols(rep))-column crosstab matrix:
         cols 1 to cols(x) give the row codes; the last cols(rep) columns give
         the frequency counts for that combination of codes and rep counts.
*/
proc crosstab(x,ncat,codes,rep,names);
local t0,k,e,s,t,i,u,mint,mult,rowcode;
t0 = date;
/* no. of vars. to be crosstabulated: */ k = cols(x);
if ncat[1]; /* use user-supplied codes: */
   if cols(x) ne rows(ncat); "INPUT ERROR IN CROSSTAB.G:";;
      " Columns of input matrix not equal to rows of ncat"; end; endif;
   if sumc(ncat) ne rows(codes); "INPUT ERROR IN CROSSTAB.G:";;
      " No. of categories not equal to no. of codes."; end; endif;
    /* ends of code strings in codes: */ e = cumsumc(ncat);
    /* starts of code strings in codes: */ s = e - ncat + 1;
    /* build columns of t giving code values: */
    i = 1; t = codes[1:e[1]];
    do until i ge k; i = i + 1;
       t = (t.*.ones(ncat[i],1))~(ones(rows(t),1).*.codes[s[i]:e[i]]);
    endo;
else; /* use all unique values as the codes (slower running): */
    /* build the columns of t that give the row values: */
    i = 1; t = unique(x[.,1],1);
    do until i ge k; i = i + 1;
       u = unique(x[.,i],1);
       t = (t.*.ones(rows(u),1))~(ones(rows(t),1).*.u);
     endo;
endif;
mint = minc(t);
/* column multipliers for summation of columns 1 to k-1
   into single code column: */
if k gt 1;
   mult = rev(cumprodc(rev(maxc(t[.,2:k])-mint[2:k]+1)))'~1;
   else; mult = 1;
endif;
/* single code for each table row: */ rowcode = sumr((t-mint').*mult);
/* build last column of t, giving count for each combination of values: */
if rows(rep) eq 1; /* no rep counts: */
   t = t~counts(sumr((x-mint').*mult),rowcode);
   else; /* use rep counts: */
   i = 0;
   do until i ge cols(rep); i = i + 1;
      t = t~countwts(sumr((x-mint').*mult),rowcode,rep[.,i]);
   endo;
endif;
if 0$+names[1];
   "Codes and frequency counts for ";$names'~"Count"; t;
   print; "Total count: ";; sumc(t[.,(cols(x)+1):cols(t)])';
          "Total rows : ";; rows(t);
   "Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(t);
endp;

