/* Test program: */ x = seqa(.1,.1,20); format /rdn 6,1;
x = x~(1+x); x~permuter(x)~permuter(x); end;
/* This proc ranomly permutes (shuffles) the rows of the matrix x. */
proc permuter(x);
retp(x[sortindc(rndu(rows(x),1)),.]);
endp;
