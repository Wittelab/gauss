/* This proc reports results from a chi-squared test.
   Inputs: x = test statistic
           df = degrees of freedom
   Output: p = P-value
*/
proc chreport(x,df);
local p; p = cdfchic(x,df);
"statistic = "$+ftos(x,"%*.*lf",3,2)$+", df = "$+ftocv(df,1,0)$+", P = ";;
format /rds 8,6;
if p ge .1; format 4,2; elseif p ge .01; format 5,3;
   elseif p ge .001; format 6,4; elseif p ge .0001; format 7,5; endif;
p; format 8,2;
retp(p);
endp;

