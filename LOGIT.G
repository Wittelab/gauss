/* This proc takes the logit of a number p between 0 and 1 */
proc logit(p);
retp( -ln(1/p - 1));
endp;
