reset
set size square
set grid
set ticslevel 0
set xr[0.01:100]
y_min = 1
y_max = 100
set yr[y_min:y_max]
set y2r[0.001:1]
set y2tics
set logscale xy
set logscale y2

set term png
dirname = './data_0_0'
#dirname = './cold_collapse'
set output dirname.'/time_radius_dens.png'

set xl 'radius[kpc]' font 'Arial,15'
set yl 'time[{/Symbol m}sec.]' font 'Arial,15'
set y2l 'density' font 'Arial,15'

gyou=1 #行番号を指定
retsu=4 #列番号を指定
loadname = dirname.'/tree_walk.dat'
xmax=system("cat " . loadname . " | awk \'NR==" . gyou ."{print $" . retsu . "}\'") # 特定の行列から値を拾う
set arrow 1 from xmax,y_min to xmax,y_max nohead lc rgb "blue" lw 1
retsu=5
t1 = system("cat " . loadname . " | awk \'NR==" . gyou ."{print $" . retsu . "}\'") # 特定の行列から値を拾う

f(x)=(3/(4*pi))/(sqrt((1+(x**2)))**5)
g(x)=real(t1)

pl dirname.'/data_0.dat' u (sqrt(($1)*($1)+($2)*($2)+($3)*($3))):($10*1000000) w dot notitle,\
f(x) w l title 'Analytical solution' axes x1y2, g(x) t ' ' lc rgb "blue"
