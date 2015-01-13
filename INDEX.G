/* Test program:  let listn = age sex "sex" state race income;
let names = sex educ race age "sex"; $listn'; $names'; index(names,listn)';
/* contrast to: */ indcv(names,listn)';
end; */
/* This proc returns the indices of the character variables in names
   as they appear in the character vector listn.
Inputs: names,listn = vectors of character variables
Output: vector of indices (positions) of names in listn
Notes: Unlike indcv in GAUSS,
       1)if names[i] does not appear in listn,
         element i of the returned vector will be zero.
       2) listn cannot have repeats
       3) the matching is case sensitive (e.g., Age ne age)
*/
proc index(names,listn);
if rows(unique(vec(listn),1)) lt rows(vec(listn));
   "INDICES.G INPUT ERROR: REPEATS IN NAME LIST."; end; endif;
retp((vec(0$+names) .eq vec(0$+listn)')*seqa(1,1,rows(vec(listn))));
endp;

