/*
/* Test program: */
 let codes = 1 2 3 7.1 8.2 9.3; let names = "blood pr" "race"; let s = 1 4;
 _field = 3; _prec = 1;
 names = names|nindc(names,codes,s);
 $names';
 end;
*/
/* This proc takes the names "names" of columns of a matrix x and returns
   a vector iname of names for the columns of the indicator matrix for
   x built using the values in codes. Companion to INDC.G.
   Inputs: names = the names of the x columns
           codes = vector of values for x columns that are to be named
           s = vector of indices at which the set of codes for column
               k of x begin.
           Ordinarily, rows(s) = rows(names). If s = 1 and rows(names) > 1,
           the same values (those in codes) will be used for each x column
   Output: iname, the indicator names made by taking up to 7 characters from
           name and appending the codes at the end.
   The codes are assumed to be be integers from 0 to 9, but this
   can be modified by specifying different values for these globals:
   _field = field width for the appended codes
   _prec = precision for appended codes
*/proc nindc(names,codes,s);
/* default field width and precision for appended codes: */
declare matrix _field = 1;
declare matrix _prec = 0;
local k,e,inames;
if rows(s) eq 1 and rows(names) gt 1; /* expand s and codes: */
   s = seqa(1,rows(codes),rows(names));
   codes = ones(rows(names),1).*.codes;
endif;
/* define endpoint(s) e: */
if rows(s) gt 1; e = (s[2:rows(s)] - 1)|rows(codes);
   else; e = rows(codes); endif;
k = 0; inames = {};
do until k ge rows(names); k = k + 1;
   inames =
inames|(0$+strsect(names[k],1,8-_field)$+ftocv(codes[s[k]:e[k]],_field,_prec));
endo;
retp(inames);
endp;
