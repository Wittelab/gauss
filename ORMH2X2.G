/* This proc does Mantel-Haenszel OR analysis of a series of 2-by-2 tables,
   provided there are no zero cells in the crude table.
   am and cm must have the same number of columns.
   If am and cm have multiple columns, it cycles through those columns,
   doing a separate analysis for each pair am[.,j], cm[.,j].
   Inputs: am = matrix of exposed cases
           bm = vector of unexposed cases
           cm = matrix of exposed noncases
           dm = col vector of unexposed noncases
           names = col vector of names of row, column, & stratifying variables
                   (set to zero if no output desired)
                   -- must have a name for each column of am
   Outputs: ormh = Mantel-Haenszel odds ratio
            lse = Robins-Breslow-Greenland SE for ln(ormh) (Biometrics,1986)
            pval = J df global Mantel-Haenszel chi-squared P-value
         -- NOTE: If there is a zero cell in the crude table, (-1,-1,p) is
                  returned, where p is the null hypergeometric probability of
                  the stratified data (same as two-sided mid-p value)
   Global: _clevel = confidence level for limits (default = 0.95)
*/
proc (3) = ormh2x2(am,bm,cm,dm,names);
local cam,j,ormh,lse,pval,a,b,c,d,ctab,
      m1,m0,n1,n0,n,stab,e,vpt,vpc,chi,r,s,rs,ss,p,q,k;
declare matrix _clevel = .95;
format /rds 8,4; names = 0$+names; cam = cols(am);
if cols(cm) ne cam; "ERROR IN ORMH2x2.G: cols(am) ne cols(cm)"; end; endif;
if names[1]; format 1,0;
   "ORMH2x2.G";; if cam gt 1; ": ";; cam;;" exposure levels";; endif; print;
   if rows(names) lt 1+cam+(rows(am) gt 1);
      " -- ERROR: number of names is ";; rows(names);;
      ", should be at least ";; 2+cam; end;
   endif; format 8,2;
endif;
j = 0; ormh = -ones(1,cam); lse = zeros(1,cam); pval = ones(1,cam);
do until j ge cam; j = j + 1; a = am[.,j]; b = bm; c = cm[.,j]; d = dm;
   if (a'd) eq 0 and (b'c) eq 0;
      if names[1];
         (""$+names[1+j]$+": No nonzero diagonal or off-diagonal products.");
      endif;
      continue;
   endif;
   /* crude table: */ ctab = sumc(a~b~c~d);
      if not prodc(ctab);
         p = prodc(phypnull(a~b~c~d));
         if names[1]; "ORMH2x2.g: The crude table for "$+names[1+j]$+" is ";;
            format 6,0; ctab';
            "The null hypergeometric probability of the stratified data is ";;
            format 8,5; p; "No large-sample statistics are appropriate.";
         endif;
         continue;
      endif;
   /* eliminate zero-margin tables: */
      stab = (a+b).*(c+d).*(a+c).*(b+d);
      stab = selif(a~b~c~d,stab .gt 0);
      a = stab[.,1]; b = stab[.,2]; c = stab[.,3]; d = stab[.,4];
   /* margins: */ m1 = a+b; m0 = c+d; n1 = a+c; n0 = b+d; n = n1+n0;
   /* odds ratio calculations: */ r = a.*d./n; s = b.*c./n;
      rs = sumc(r); ss = sumc(s);
      if ss; ormh[j] = rs/ss;
         if rs; /* OR SE calculations: */ p = (a+d)./n; q = (b+c)./n;
            lse[j] = sqrt((sumc(p.*r)/(rs^2)
                           + sumc(p.*s+q.*r)/(rs*ss) + sumc(q.*s)/(ss^2))/2);
            else; if names[1]; "WARNING: M-H odds ratio is zero."; endif;
                  continue;
         endif;
         else;
           if names[1]; "WARNING: M-H odds ratio infinite or undefined"; endif;
           continue;
      endif;
   /* null expected a-cell: */ e = m1.*n1./n; vpt = e.*n0./n;
   /* Null a-cell variances: */ vpc = vpt.*m0./(n-1);
   /* M-H (score) statistic: */
      if sumc(vpc); chi = (ctab[1] - sumc(e))./sqrt(sumc(vpc));
         else; chi = 0; endif; pval[j] = cdfnt(chi);
   if names[1]; /* print results: */
      (" -- Mantel-Haenszel analysis of "$+names[1]$+" and "$+names[1+j]);;
      if rows(names) gt 1+cam; format /rds 8,4;
         ","; "    stratified by: ";; $names[2+cam:rows(names)]';;
         "(";; format 1,0; rows(a);;"strata)";;
      endif; ".";
      if rs*ss; local z; z = cdfni(1-(1-_clevel)/2); z = 0~(-z)~z;
      "    M-H OR estimate, & ";; format /rdn 1,0; 100*_clevel;; "% limits:";;
         format /rds 9,3; ormh[j]*exp(z*lse[j]);
      endif;
      "    Mantel-Haenszel chi statistic & 2-sided p:";; format 8,2;
      chi;; ", ";; format 7,4; pval[j];
      "    crude OR and 2x2 table    :";;
      if ctab[2]*ctab[3]; format 9,3; ctab[1]*ctab[4]/(ctab[2]*ctab[3]);;
         else; "INFINITE ";; endif;
      format 6,0; "   ";; ctab'; format 9,3;
   endif;
endo;
if cam gt 1 and prodall(ormh .gt 0); /* compute global P-value: */
   m1 = sumr(am~bm); m0 = sumr(cm~dm); n1 = am+cm; n0 = bm+dm;
   n = m1+m0; n = n + (n .eq 0); e = m1.*n1./n;
   vpt = (n1./n).*m0./(n-(n .gt 1)); vpc = zeros(cam,cam); k = 0;
   do until k ge rows(n); k = k + 1; vpc = vpc - e[k,.].*vpt[k,.]'; endo;
   vpc = diagrv(vpc,sumc(e.*m0./(n-(n .gt 1)) - vpt)); rs = sumc(am-e);
   trap 0; pval = cdfchic(rs'invpd(vpc)*rs,cam);
   if scalerr(pval); pval = -1; trap 1;
      "ORMH2x2: GLOBAL COVARIANCE MATRIX NOT POSITIVE DEFINITE";
      elseif names[1]; format 7,4;
      (""$+ftocv(cam,1,0)$+" df global Mantel-Haenszel P-value = ");; pval;
   endif;
endif;
retp(ormh,lse,pval);
endp;

