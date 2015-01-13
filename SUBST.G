/* Test program: format /rdn 7,0; let v = 7 8 9;
let s[3,2] = 1 2  2 3  3 4; subst(-ones(3,4),s[.,.],v); end; */
/* This proc substitutes the value from v
   into the matrix x in the columns given in s. */
proc subst(x,s,v);
local j,sx;
j = 0; sx = x;
do until j ge cols(s); j = j + 1;
   sx = substute(sx,s[.,j] .eq seqa(1,1,cols(x))',v[.,min(j,cols(v))]);
endo;
retp(sx);
endp;
