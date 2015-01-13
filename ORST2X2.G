/* This proc does standardized OR analysis of a series of 2-by-2 tables,
   using the total population distribution (as imputed from the cases)
   as standard, provided there are no zero cells in the strata.
   Not recommended if there are many small cells within strata.
   Inputs: a = column vector of exposed cases
           b = col vector of unexposed cases
           c = col vector of exposed noncases
           d = col vector of unexposed noncases
           names = col vector of names of row, column, & stratifying variables
                   (set to zero if no output desired)
   Outputs: orst = odds ratio standardized to total
            lse = SE for ln(orst)
            1 df Mantel-Haenszel chi-squared p-value
         -- NOTE: If there is a zero cell in the crude table, (-1,-1,p) is
                  returned, where p is the null hypergeometric probability of
                  the stratified data (same as two-sided mid-p value)
   Global: _clevel = confidence level for limits (default = 0.95)
*/
proc (3) = orst2x2(a,b,c,d,names);
local ctab,m1,m0,n1,n0,n,e,vpt,vpc,chi,r,s,rs,ss,orst,p,q,lse;
declare matrix _clevel = .95; format /rds 8,4;
names = 0$+names;
if sumc((a.*b.*c.*d) .eq 0);
   if names[1]; "ORST2X2.G: Zero cells were found"; endif; retp(-1,-1,1);
endif;
/* crude table: */ ctab = sumc(a~b~c~d);
/* margins: */ m1 = a+b; m0 = c+d; n1 = a+c; n0 = b+d; n = n1+n0;
/* odds ratio calculations: */ r = a.*m0./c; s = b.*m0./d;
   rs = sumc(r); ss = sumc(s);
   if ss; orst = rs/ss;
      if rs; /* OR SE calculations: */ lse =
         sqrt(((r^2)'(1/a+1/c+1/m0))/(rs^2) + ((s^2)'(1/b+1/d+1/m0))/(ss^2));
         else; lse = 0;
           if names[1]; "WARNING: Standardized odds ratio is zero"; endif;
      endif;
      else; orst = -1; lse = 0;
           if names[1];
              "WARNING: Standardized odds ratio infinite or undefined"; endif;
   endif;
/* null expected a-cell: */ e = m1.*n1./n; vpt = e.*n0./n;
/* Null a-cell variances: */ vpc = vpt.*m0./(n-1);
/* M-H (score) statistic: */
   if sumc(vpc); chi = (ctab[1] - sumc(e))./sqrt(sumc(vpc));
      else; chi = 0; endif;
if names[1]; /* print results: */ print;
   ("ORST2X2.G: Standardized analysis of "$+names[1]$+" and "$+names[2]);;
   if rows(names) ge 3; format /rds 8,4;
      ","; " stratified by: ";; $names[3:rows(names)]';;
      "(";; format 1,0; rows(a);;"strata)";;
   endif; ".";
   if rs*ss; local z; z = cdfni(1-(1-_clevel)/2); z = 0~(-z)~z;
      "    Std OR estimate, & ";; format /rdn 1,0; 100*_clevel;; "% limits:";;
      format /rds 9,3; orst*exp(z*lse);
   endif;
   "    Mantel-Haenszel chi statistic & 2-sided p:";; format 8,2;
      chi;; if chi; format 8,5; cdfnt(chi); else; ", 1.0"; endif;
   "    crude OR and 2x2 table    :";;
   if ctab[2]*ctab[3]; format 9,3; ctab[1]*ctab[4]/(ctab[2]*ctab[3]);;
      else; "INFINITE ";; endif;
   format 6,0; "   ";; ctab'; format 9,3;
endif;
retp(orst,lse,cdfnt(chi));
endp;
