proc (0) = descm(x);
local m;
" min, median, mean, max:";; m = minc(x)~medianc(x)~meanc(x)~maxc(x); 
if cols(x) .gt 1; print; endif; m;
endp;