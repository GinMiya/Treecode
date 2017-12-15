reset
set size square
set view equal xy
set grid
set ticslevel 0
set xr[-1:1]
set yr[-1:1]
#set logscale z

set term png
#unset ztics

dirname = './data_0_0'
set output dirname.'/time_radius_2d.png'

set xl 'x[kpc]' font 'Arial,15'
set yl 'y[kpc]' font 'Arial,15'
set yl rotate by 90
set zl 'time[{/Symbol m}sec.]' font 'Arial,15'

gyou=1 #行番号を指定
retsu=1 #列番号を指定
loadname = dirname.'/tree_walk.dat'
x1=system("cat " . loadname . " | awk \'NR==" . gyou ."{print $" . retsu . "}\'") # 特定の行列から値を拾う
retsu=2
y1=system("cat " . loadname . " | awk \'NR==" . gyou ."{print $" . retsu . "}\'") # 特定の行列から値を拾う
retsu=5
t1=system("cat " . loadname . " | awk \'NR==" . gyou ."{print $" . retsu . "}\'") # 特定の行列から値を拾う
set label 1 point pt 7 ps 1 lc rgb "red" at x1,y1

#set view 0,0,,
spl dirname.'/data_0.dat' u 1:2:10 w dot notitle
