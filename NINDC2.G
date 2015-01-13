/* This proc takes the names "names" of columns of a matrix x and returns
   a vector iname of names for the columns of the indicator matrix for
   x, naming them names1,names2,etc., where names gets truncated at the
   seventh character. Companion to INDC2.G.
   Inputs: nind = vector of number of indicators for each column
           names = the names of the x columns
   Output: iname, the indicator names made by taking up to 7 characters from
           name and appending the codes at the end.
   The numberof indicators for each column is assumed to be less than 9,
   but this can be modified by specifying different values for these globals:
   _field = field width for the appended codes
   _prec = precision for appended codes
*/
/* Test program:
let nind = 3 4; let names = "blood pr" "race";
names = names|nindc2(nind,names); $names'; end; */
proc nindc2(nind,names);
/* default field width and precision for appended codes: */
declare matrix _field = 1;
declare matrix _prec = 0;
local k,e,inames;
k = 0; inames = {};
do until k ge rows(names); k = k + 1;
   inames = inames|
(0$+strsect(names[k],1,8-_field)$+ftocv(seqa(1,1,nind[k]),_field,_prec));
endo;
retp(inames);
endp;
