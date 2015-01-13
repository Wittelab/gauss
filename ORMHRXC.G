/* This proc does Mantel-Haenszel and exact
   odds-ratio analysis of a series of r-by-c tables,
   treating the first row and column as the reference levels
   (reference levels may be changed by setting globals; see below)
   Inputs: nr = no. of rows of each table
           t = stacked r-by-c tables (counts)
           rnames = col vector of row names
                   (set to zero if no output desired)
           cnames = col vector of col names
                   (set to zero if no output desired)
           snames = col vector of stratifying variable names
                   (set to zero if there is only one stratum)
   Outputs: ormh = r-1 by c-1 table of Mantel-Haenszel odds ratios
            sel = table of Robins-Breslow-Greenland SEs for the ln(ormh)
            pval = table of 1 df M-H chi-squared p-values.
   Globals: _rf = reference row (default is 1)
            _cf = reference column (default is 1)
            _clevel = confidence level for limits (default = 0.95)
            _mincc = minimum crude cell size for use of M-H statistics
                     (default = 5)
            _minp = minimum M-H p-value for use of M-H stats (default = 0.05)
   NOTE 1: If any cell in the crude table is less than _mincc or
           the M-H p is less than _minp, exact mid-p statistics
           will be returned for that cell instead of M-H statistics.
           The minimum crude-cell criterion is not sufficient, however (see
           Mantel & Fleiss, AJE 1981), and can be increased by raising _mincc.
           If both _mincc & _minp are set to 0, no exact stats will be used.
   NOTE 2: If there is a zero cell in the crude table, (-1,-1,p) is
           returned for that cell, where p is the null hypergeometric
           probability of the corresponding stratified 2-by-2 data,
           equal to the two-tailed (by doubling rule) mid-p value.
*/
proc (3) = ormhrxc(nr,t,rnames,cnames,snames);
declare matrix _rf = 1; declare matrix _cf = 1; declare matrix _clevel = .95;
declare matrix _mincc = 5; declare matrix _minp = .05;
local r0,rf,nc,r,c,ormh,sel,pval,orex,orl,oru,fl,fu,z,
      blank,infin,undef,zed,v;
if fmod(rows(t),nr);
   "ORMHRxC.G: NUMBER ROWS SPECIFIED INCONSISTENT WITH INPUT TABLE"; end;
endif;
r0 = seqa(0,nr,rows(t)/nr); rf = r0 + _rf; nc = cols(t);
ormh = -ones(nr,nc); ormh[_rf,.] = ones(1,nc); ormh[.,_cf] = ones(nr,1);
pval = ormh; sel = ormh; sel[_rf,.] = zeros(1,nc); sel[.,_cf] = zeros(nr,1);
orex = ormh; orl = ormh; oru = ormh; r = r0;
do until r[1] ge nr; r = r + 1;
  if r[1] ne _rf; c = 0;
    do until c ge nc; c = c + 1;
      if c ne _cf;
         { ormh[r[1],c],sel[r[1],c],pval[r[1],c] } =
           ormh2x2(t[r,c],t[rf,c],t[r,_cf],t[rf,_cf],0);
         if (pval[r[1],c] lt _minp and pval[r[1],c] ge 0)
            or sumall((t[r,c _cf]~t[rf,c _cf]) .lt _mincc);
            { orex[r[1],c],orl[r[1],c],oru[r[1],c],fl,fu,pval[r[1],c] } =
              exact2x2(t[r,c],t[rf,c],t[r,_cf],t[rf,_cf],0);
         endif;
      endif;
    endo;
  endif;
endo;
if 0$+rnames[1]; /* print results: */
   print; "Proc ORMHRxC.G: Tabular analysis of the following data,";
   if 0$+snames[1]; "   stratified on ";;
     if rows(snames) gt 1; format /rdc 9,0; $snames[1:rows(snames)-1]';; endif;
     format /rds 9,0; $snames[rows(snames)];
   endif;
   "  Crude counts, Mantel-Haenszel odds-ratio estimates,";
   "  Robins et al. ";; format /rdn 3,0; 100*_clevel;; "% limits, and 2-sided";
   "  Mantel-Haenszel p-values, or exact mid-p estimates and 2-sided";
   "  (twice the smallest one-sided) mid-p values where noted by *.";
   "  A -1 will be printed where an exact limit failed to converge.";
   z = cdfni(1-(1-_clevel)/2); z = 0~(-z)~z;
   let blank = " "; let infin = "infinite";
   let undef = "     not" "defined "; let zed = "0";
   c = 0;
   do until c ge nc; c = c + 1;
     if c ne _cf; r = 0;
      print; "  FOR COLUMN "$+cnames[c]$+":"; format /rds 8,0;
"             ";;$cnames[c _cf]';;"  odds ratio    lower & upper CL    p-val";
      do until r ge nr; r = r + 1;
         "    ";; format /lds 8,0; $rnames[r];;
         format /rds 8,0; sumc(t[r0+r,c _cf])';; format 10,2;
         if r ne _rf;
            if orex[r,c] gt 0; orex[r,c]~orl[r,c]~oru[r,c];;"* ";;
               elseif ormh[r,c] gt 0; ormh[r,c].*exp(z*sel[r,c]);;"  ";;
               elseif orl[r,c] gt 0; $infin;; orl[r,c];; $infin;;"* ";;
               elseif oru[r,c] gt 0; $zed[1 1];; oru[r,c];;"* ";;
               else; $(undef'~blank);;"* ";;
            endif;
            v = pval[r,c]; format /ldn 6,4;
            if v ge .1; format 6,2; elseif v ge .01; format 6,3; endif; v;
            else; "  (reference row)";
         endif;
      endo;
     endif;
   endo;
endif;
retp(ormh,sel,pval);
endp;
