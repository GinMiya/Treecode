if(exist ("n") == 0 || n<0) n = n0

title(n) = sprintf("t = %.2f",0.01*n)
set title title(n) font 'Arial,20'
loadname = dirname.sprintf("data_%d.dat",n)
pl loadname u 7:8 w dot notitle

n = n+dn
if(n<nm)reread
