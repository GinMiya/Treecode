reset
set ticslevel 0
set size square
set grid
set logscale x
#set xr[10**(-9):10**(0)]
set xr[10**(-7):10**(3)]
y_min = 0
y_max = 100
set yr[y_min:y_max]
set key left
set term png enhanced

dirname = "./data_0_0"
#dirname = "./cold_collapse"
filename = dirname.'/acc_error_percentile.png'
loadname = dirname.'/acc_error.dat'


set output filename
set format x "10^{%L}"

set xl 'Force error magnitude' font "Arial,15"
set yl 'Percentile' font 'Arial,15'

gyou  = 1 #行番号を指定
retsu = 1 #列番号を指定
xmax=system("cat " . loadname . " | awk \'NR==" . gyou ."{print $" . retsu . "}\'") # 特定の行列から値を拾う
#set label at graph xxmax,0.5 sprintf("max=%f",xmax) font 'Arial,10'
#set arrow 1 from xmax,y_min to xmax,y_max nohead lc rgb "forest-green" lw 1

f(x)=99.0
g(x)=90
# plot loadname index 1 u 2:1 w l notitle,f(x) w l title '99%',g(x) w l title '90%'
plot loadname u 2:1 w l title 'Warren{/text \46}Salmon', f(x) w l title '99%',g(x) w l title '90%'\
#'./data_0_0/acc_error_multi.dat' u 2:1 w l title 'BH Tree',\
