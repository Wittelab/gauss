/* This proc prints a numeric matrix with beginning row names.
NOTE: Format controlled by calling program.
Inputs: names = row names, x = numeric matrix
To yield flush left margin, pad out names to 8 characters with trailing blanks.
*/
proc (0) = nameprin(names,x);
local i; i = 0; do until i ge rows(x); i = i + 1; $names[i];; x[i,.]; endo;
retp; endp;
