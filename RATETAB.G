/* This proc crosstabulates the ratios sumc(y)'/sumc(d) within subcategories of
   the variables in x.
   Inputs: x = crosstab (stratification) variables
           y = numerator (case) indicators or counts (may have multiple cols)
           d = person-time for person-time rates, total count for risks,
               noncases for odds or case-control ratios;
               if a scalar, d will be redefined as a vector times that scalar.
           ncat = vector giving number of code values to be used for each
                  column in x -- if 0, program will use all unique values
                  in a column as the column codes.
           codes = vector of length sumc(ncat) giving in succession the code
                   values for each variable in x (will be ignored if ncat = 0)
           rep = repetition count (0 if none)
           xnames, ynames, dname, rname = names of x cols, y cols, d, & y/d
                                   (set to 0 for no printed output)
   Outputs: t = prodc(ncat)-row by (cols(x)+2*cols(y)+1)-column crosstab matrix:
            cols 1 to cols(x) give the row (stratum) codes; the other columns
            give sumc(y)', sumc(d), and sumc(y)'/sumc(d) for
            that stratum (combination of codes).
   Globals: _cl = if this is nonzero, y has only one column, and names are
                  supplied, print 95% confidence limits for the ratios
                  (default 0); 1 for rates, 2 for risks, 3 for odds.
            _dzero = code for zero divide; default -1.
            _delzero = 1 to delete rows with zero denominators from output
                       and return; default 0.
   OUTPUT NOTES:
         1) When d is zero, a "zero-divide" value _dzero wiil be printed for
            the ratio. This global variable has default value -1.
         2) You may adjust the formats of the header and table body by setting
            the printfm matrices _fmth and _fmt in the calling program.
*/
proc ratetab(x,y,d,ncat,codes,rep,xnames,ynames,dname,rname);
declare matrix _cl = 0;
declare matrix _dzero = -1; declare matrix _delzero = 0;
declare matrix _fmt = 0; declare matrix _fmth = 0;
local t0,nx,ny,t,r,z,und,ind,v,lims;
t0 = date; nx = cols(x); ny = cols(y);
if rows(d) eq 1; d = d*ones(rows(x),1); endif;
if rows(rep) eq 1; /* no rep count: */
   t = crosstab(x,ncat,codes,y~d,0);
   else; /* use rep count: */
   t = crosstab(x,ncat,codes,rep.*(y~d),0);
endif;
y = t[.,(nx+1):(nx+ny)]; d = t[.,(nx+1+ny)]; und = d .eq 0;
/* ratios: */ r = (1-und).*y./(d + und) + und.*_dzero;
ind = seqa(1,1,rows(t));
if _delzero; /* index nonzero-denominator rows: */
   ind = delif(ind,und); endif;
if 0$+xnames[1]; local fmt,codefmt,fmth,ratefmt; print;
   " Proc ratetab.g: Tabulation of stratum-specific counts & "$+rname$+".";
   if _fmt[1,1] eq 0; /* build default output formats: */
      let codefmt[1,3] = "s " 8 8;
      codefmt = ones(nx+ny,1).*.codefmt;
      let fmth[2,3] = "s " 8 8 "s " 9 9;
      fmth = codefmt|fmth[1,.]|(ones(ny,1).*.fmth[2,.]);
      let codefmt[1,3] = "lf " 8 0;
      codefmt = ones(nx+ny,1).*.codefmt;
      let fmt[1,3] = "lf " 8 0; let ratefmt[1,3] = "lf " 9 2;
      ratefmt = ones(ny,1).*.ratefmt;
      fmt = codefmt|fmt|ratefmt;
      if ny eq 1 and _cl; if _cl eq 1; fmt[nx+ny+1,3] = 2; endif;
         let codefmt[2,3] = "s " 9 9 "s " 8 8; fmth = fmth|codefmt;
         let ratefmt[2,3] = "lf," 9 2 "lf " 8 2; fmt = fmt|ratefmt;
      endif;
      else; fmt = _fmt; fmth = _fmth;
   endif;
   rname = ones(ny,1)'rname;
   if ny eq 1 and _cl;
      z = (-1.96)~1.96; und = y .le 0; v = 1/(y+und);
     if _cl eq 2; und = und .or d .le 0; v = v - 1/(d+und);
        elseif _cl eq 3; und = und .or (d-y) .le 0; v = v + 1/(d-y+und);
           endif;
      lims = r.*exp(z.*sqrt(v)); if _cl eq 3; lims = lims./(1+lims); endif;
      r = r~((1-und).*lims + _dzero*und);
      rname = rname~"95% ll"~"95% ul";
   endif;
   call printfm(xnames'~(ynames')~dname~rname,zeros(1,rows(fmth)),fmth);
   print; call printfm(t[ind,.]~r[ind,.],ones(1,rows(fmt)),fmt);
   " Totals & crude rates: "; v = sumc(d) eq 0;
   call printfm((" ".*.ones(1,nx))~(sumc(y~d)')~
           ((1-v)*sumc(y)'/(sumc(d)+v) + _dzero*v),
           zeros(1,nx)~ones(1,2*ny+1),
           fmth[1:nx,.]|fmt[nx+1:nx+2*ny+1,.]);
   print;
   if sumc(und) and not _delzero;
      format 4,0; "A value of ";;_dzero;;" denotes a zero denominator."; endif;
   " Total run time: ";; etstr(ethsec(t0,date));
endif;
retp(t[ind,.]~r[ind,.]);
endp;
