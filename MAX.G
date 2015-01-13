/* This proc returns the elementwise maximum of x and y */
proc max(x,y);
retp((x .gt y).*x + (x .le y).*y);
endp;
