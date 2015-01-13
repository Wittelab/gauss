/* Test program:
let v = 1 2  1 1 1 4 5 5 6;
format /rdn 4,0; v'; { ic,lnt } = changept(v); ic'; lnt';
end;
*/
/* This proc locates the change points in a sequence (vector) v.
   Input: v = vector of numbers.
   Outputs: ic = indices of change points in v,
            length of each run of constant value in v.
*/
proc (2) = changept(v);
local nr,ic;
nr = rows(v);
/* ic = v shifted down 1 index: */ ic = (v[1]-1)|v[1:(nr-1)];
/* ic = indices of change points: */ ic = selif(seqa(1,1,nr),ic .ne v);
retp(ic,(ic[2:rows(ic)]|(nr+1)) - ic);
endp;

