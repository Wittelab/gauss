/* This proc creates names for interaction (product) terms for regressions
   for two separate sets of regressors cross-multiplied, i.e., for x1*~x2.
   Inputs: names1, names2 = names of x1,x2 columns
   Outputs: namint = name vector for products (e.g., "name1*name2")
   WARNING: Because names can have only 8 characters, your input names will be
            truncated to 3 or 4 characters - be sure you can distinguish your
            variables based on only 3 characters.
*/
proc ninter2(names1,names2);
local k; names1 = 0$+names1; names2 = 0$+names2;
k = 0;
do until k ge rows(names1); k = k + 1;
   names1[k] = (0$+strsect(names1[k],1,4));
endo;
retp((names1.*.ones(rows(names2),1))$+"*"$+(ones(rows(names1),1).*.names2));
endp;
