/* test program: negpower(-(2|8|4),0|(1/3)|(-1/2)); */
/* This proc returns the simple spline sign(x).*abs(x)^p,
   which is useful for representing fractional powers
   of potentially negative variables in dose-response functions. */
proc negpower(x,p);
retp(sign(x).*((abs(x)+(x .eq 0))^p));
endp;

