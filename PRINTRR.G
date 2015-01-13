/* This proc prints a matrix of estimated ratios and the corresponding matrices
   of interval limits to 2-digit accuracy, in one of two formats:
   estimate and interval side-by-side, or estimate and interval above/below.
Inputs: rnames = column vector of row names for r (0 if no names to print)
        cnames = column vector of column names for r (0 if no names to print)
        r, rl, ru = ratios, lower limits, upper limits
        ab = 0 for side-by-side, 1 for above/below
NOTE: If an element of r, rl, or ru is negative,
      a blank will be printed in its place
*/
proc (0) = printrr(rnames,cnames,r,rl,ru,ab);
local nc,bl,j,i; cnames = 0$+cnames; rnames = 0$+rnames; format /rdn 8,3;
nc = cols(r); bl = "        ";
if cnames[1]; bl;; "  ";; if ab; "  ";; endif; j = 0;
   do until j ge nc; j = j + 1;
      $cnames[j];; bl;; if not ab and j lt nc; bl;; endif;
   endo; print;
endif;
i = 0;
do until i ge rows(r); i = i + 1; format 8,2;
   if rnames[1]; $rnames[i];; endif; if ab; "   ";; endif;
   j = 0;
   do until j ge nc; j = j + 1; format 8,2;
      if r[i,j] ge 0;
         if ab; r[i,j];; bl;;
            else;  format 7,2; r[i,j];;
            " (";; format 5,2; rl[i,j];; ",";; format 7,2; ru[i,j];;")";;
            if j lt nc; " ";; endif;
         endif;
         elseif j lt nc; bl;; bl;; if not ab; bl;; endif;
      endif;
   endo; print;
   if ab; if rnames[1]; bl;; endif;
      "  ";; j = 0;
      do until j ge nc; j = j + 1;
         if r[i,j] ge 0;
            "(";; format 5,2; rl[i,j];; ",";; format 7,2; ru[i,j];;") ";;
            else; bl;; bl;;
         endif;
      endo; print;
   endif;
endo; print;
retp;
endp;
