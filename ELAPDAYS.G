/* This proc takes contiguous-element index date(s) and returns
   days elapsed from reference date(s) to the index date(s).
   Inputs: d = date(s) to be converted
           r = reference date or date vector (in same format as d)
           order = 1 if month first, day second
                   0 if day comes first, month second
   NOTES: 1) d MUST be in 6-digit format with a leading zero on the middle
          component when that component is 1-9, e.g., if month is in middle,
          months must be coded 01, 02, etc.
          2) years are assumed to have two digits and all be within the same
          century.
          3) No index date may be earlier than the corresponding ref. date.
   Output: t = elapsed time(s), in days.
*/
/* Test program:
let d1 = 042092 032093 032293; r = 032192; order = 1;
dt1 = elapdays(d1,r,order); d1~dt1;
let d2 = 200492 200393 220393; r = 210392; order = 0;
dt2 = elapdays(d2,r,order); format /rdn 10,0; d2~dt2;
let d3 = 032293; order = 1;
dt3 = elapdays(d3,d1,order); format /rdn 10,0; d1~dt3;
end; */
proc elapdays(d,r,order);
local day,month,year,rday,rmonth,ryear,i,t,dr,de;
if rows(r) eq 1; r = r*ones(rows(d),1);
   elseif rows(d) eq 1; d = d*ones(rows(r),1);
endif;
year = d%100; ryear = r%100;
if order;
   day = floor((d%10000)/100); month = floor(d/10000);
   rday = floor((r%10000)/100); rmonth = floor(r/10000);
else;
   month = floor((d%10000)/100); day = floor(d/10000);
   rmonth = floor((r%10000)/100); rday = floor(r/10000);
endif;
i = 0; t = zeros(rows(d),1);
do until i ge rows(d); i = i + 1;
    dr = ryear[i]|rmonth[i]|rday[i];
    de = year[i]|month[i]|day[i];
    t[i] = etdays(dr,de);
endo;
retp(t);
endp;

