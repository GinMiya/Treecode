reset
set size square
set grid
set ticslevel 0
y_min = 0
y_max = 100
set yr[y_min:y_max]

set term png
dirname = './data_0_0'
#dirname = './cold_collapse'
set output dirname.'/time_percentile.png'

#set label at graph xxmax,0.5 sprintf("max=%.1f",xmax) font 'Arial,10'
gyou=1 #行番号を指定
retsu=2 #列番号を指定
loadname= dirname.'/acc_error.dat' # ファイル名を指定
xmax = system("cat " . loadname . " | awk \'NR==" . gyou ."{print $" . retsu . "}\'") # 特定の行列から値を拾う
set arrow 1 from xmax,y_min to xmax,y_max nohead lc rgb "forest-green" lw 1

set xl 'Time[μsec.]' font 'Arial,15'
set yl 'Percentile' font 'Arial,15'

pl loadname index 1 u ($3*1000000):1 smooth cspline w l notitle
