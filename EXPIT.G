/* This proc does the logistic transform on x */
proc expit(x);
retp(1/(1 + exp(-x)));
endp;

