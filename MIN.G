/* This proc returns the elementwise minimum of x and y */
proc min(x,y);
retp((x .lt y).*x + (x .ge y).*y);
endp;
