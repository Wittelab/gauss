/* This proc generates a fully crossed design matrix of variables where
   nl = vector containing number of levels of each factor in sequence
*/
proc crossmat(nl);
local nf,x,i; nf = rows(nl); x = zeros(prodc(nl),nf);
x[.,1] = seqa(0,1,nl[1]).*.ones(prodc(nl[2:nf]),1);
i = 1;
do until i ge nf-1; i = i + 1;
   x[.,i] =
   ones(prodc(nl[1:i-1]),1).*.seqa(0,1,nl[i]).*.ones(prodc(nl[i+1:nf]),1);
endo;
x[.,nf] = ones(prodc(nl[1:nf-1]),1).*.seqa(0,1,nl[nf]);
retp(x);
endp;

