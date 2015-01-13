/* This proc creates names for interaction (product) terms for regressions.
   Companion to INTER.G.
   Inputs: names = names of x columns
   Outputs: namint = name vector for products (e.g., "name1*name2")
   NOTE: If input vector has only one element the proc will return {}.
   WARNING: Because names can have only 8 characters, your input names will be
            truncated to 3 or 4 characters - be sure you can distinguish your
            variables based on only 3 characters.
*/
/* See INTER.G for test program */
proc ninter(names);
local namint,k,i;
namint = {};
if rows(names) eq 1; retp(namint); endif;
k=0;
do until k ge rows(names); k = k + 1;
   i = k;
   do until i ge rows(names); i = i + 1;
   namint =
     namint|((0$+strsect(names[k],1,4))$+"*"$+(0$+strsect(names[i],1,3)));
   endo;
endo;
retp(namint);
endp;


