/* This proc takes a bootstrap sample from the data matrix x using the proc
   rndsrep (random sampling with replacement). */
proc bootst(x);
retp(x[rndsrep(rows(x),rows(x)),.]);
endp;
