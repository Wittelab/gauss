/* Computes summary statistics for each ROW of a matrix b 
   of loglinear/logistic coeffcients from a simulation: 
*/
proc (0) = simsum(b);
local zcut,nc,pct,j,bj;
zcut = 0~(-1.96)~1.96;
print; format 8,2; nc = cols(b); 
"Antilogs of means,-/+1.96sd, medians, 2.5 and 97.5 %tiles of "$+ftocv(nc,1,0)$+" total:";
pct = ceil(nc/2)|floor(nc/40)|ceil(39*nc/40);
j = 0;
do until j ge rows(b); j = j+1;
   bj = sortc(b[j,.]',1);
   exp(meanc(bj)+zcut*stdc(bj));;"    ";; exp(bj[pct])';
endo;
"-- actual percentiles:";; 100*pct'/nc;
endp;
