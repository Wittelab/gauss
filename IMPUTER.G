/* Imputes missing values from complete records in dat via OLS.
   Inputs: dat = data matrix
           misind = missing indicator matrix for dat
           names = names of dat columns (0 for no printed output)
   Output: dat = imputed version of dat
   Globals: _const = 1 for constant (default), 0 for none
*/
proc imputer(dat,misind,names);
declare matrix _const = 1;
local nmis,x,nc,sr,sc,k,s,u,old,b;
names = 0$+names;
nmis = sumc(misind);
nc = cols(dat);
x = delif(dat,sumr(misind));
if _const; x = x~ones(rows(x),1); endif;
if names[1]; 
   "Imputing missings by OLS; complete records = "$+ftocv(rows(x),1,0);
endif;
sr = seqa(1,1,rows(dat)); 
sc = seqa(1,1,cols(x)); 
k=0;
do until k ge rows(nmis); k=k+1; 
   if not nmis[k]; continue; endif;
   s = selif(sr,misind[.,k]); 
   u = delif(sc, sc .eq k); 
   old = dat[s,k]; 
   b = invpd(x[.,u]'x[.,u])*(x[.,u]'x[.,k]);
   dat[s,k] = dat[s,u[1:nc-1]]*b[1:nc-1];
   if _const; dat[s,k] = dat[s,k] + b[nc]; endif;
   if names[1];
      ""$+names[k]$+" has "$+ftocv(nmis[k],1,0)$+" imputations:";
      old';; "becomes";;
      dat[s,k]'; print;
   endif;
endo;
retp(dat);
endp;