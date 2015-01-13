/* assigns the global names in names to the columns of input matrix x */
proc (0) = varputs(x,names);
local i,s;
if rows(names) ne cols(x); 
   "ERROR IN VARPUTS.G: rows(names) ne cols(x)"; end;
endif;
i=0;
do while i < rows(names); i = i+1;
   s=varput(x[.,i],names[i]);
endo;
endp;