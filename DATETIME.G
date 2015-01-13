/* This proc prints the date and time in the form "day mo, yr, at hr:min" */
proc (0) = datetime;
local month,d,t,ampm;
let month = "Jan." " Feb." "Mar." "Apr." "May." "June"
            "July" "Aug." "Sept" "Oct." "Nov." "Dec.";
let ampm = "am" "pm";
d = date; t = time;
format /rdn 2,0; d[3];; format 4,0; " ";; $month[d[2]];;" ";; d[1];;
" at ";; format 2,0; (t[1] % 12);; ":";; t[2];; " ";;
$ampm[((t[1] ge 12)+1)];
retp;
endp;

