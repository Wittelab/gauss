/* This procedure computes the Pearson chi-squared for a 2-way table.
   Inputs: t = 2-way table of counts
           rnames = names for the rows (0 for no printed output)
           cnames = names for the columns (0 for no printed output)
   Output: x2p = Pearson chisquared for t
           df = degrees of freedom for x2p
           p-value for x2p
*/
proc (3) = pearson(t,rnames,cnames);
local tr,e,x2p,df,p; tr = sumr(t);
/* expecteds: */ e = tr.*(sumc(t)'/sumc(tr));
x2p = sumall(((t-e)^2)./e);
df = (rows(t)-1)*(cols(t)-1); p = cdfchic(x2p,df);
if rnames[1] and cnames[1]; /* print table and results: */ local i,tot;
   i = 0; let tot = "totals";
   "For the table"; "          ";; format /rds 9,0; $cnames';;$tot;
   i = 0;
   do until i ge rows(t); i = i + 1;
      format /lds; $rnames[i];; format /rds; t[i,.];; tr[i];
   endo;
   format /lds; $tot;; format /rds; sumc(t~tr)';
   "Pearson chi-squared: ";; format 12,2; x2p;;
   " on ";; format 1,0; df;; " degrees of freedom, p =";;
   format 8,5; p; format 9,3;
endif;
retp(x2p,df,cdfchic(x2p,df));
endp;

