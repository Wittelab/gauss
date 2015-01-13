/* Builds a complete prodc(n) by rows(n) design matrix from quantitative codes:
    codes = vector of concatenated codes from slowest 
                  to fastest changing variables; rows must equal prodc(n)
    n = number of codes for each variable from slowest 
           to fastest changing variable; product must equal rows(codes)
*/
proc designq(codes,n);
local r,c,d,i,s;
r = rows(n); if r eq 1; retp(codes); endif;
c = cumprodc(n);
d = zeros(c[r],r);
d[.,1] = codes[1:n[1]].*.ones(prodc(n[2:r]),1);
i=1; s=n[1];
do until i ge r-1; i = i+1;
     d[.,i] = ones(c[i-1],1).*.codes[s+1:s+n[i]].*.ones(prodc(n[i+1:r]),1);
     s = s+n[i];
endo;
d[.,r] = ones(c[r-1],1).*.codes[s+1:s+n[r]];
retp(d);
endp;