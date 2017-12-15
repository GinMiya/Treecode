reset
set size square
set grid
set ticslevel 0
set xr[0.01:100]
y_min = 1
y_max = 1000
set yr[y_min:y_max]
set rmargin 50
#set view equal xy
set logscale xy
set logscale cb 2
set format cb "2^{%L}"

hini = 0.001
vini = 0.001
hsize=1.1
vsize=1.1
set bmargin -10

set term png
dirname = './data_0_0'
#dirname = './cold_collapse'
loadname = dirname.'/tree_walk.dat'
set output dirname.'/nbody_radius.png'

set xl 'radius[kpc]' font 'Arial,15'
#set yl 'N particle in CELLs' font 'Arial,15' rotate by 90
set origin hini,vini
set si hsize,vsize
set label 1 'N particle in CELLs' at graph 1.15, graph 0.1 font 'Arial,15' rotate by 90
set label 2 'CELL length' at screen 0.2, screen 0.42 font 'Arial,15' rotate by 90
set colorbox user origin 0.1,0.32 size 0.03,0.4
unset ztics

gyou=1 #行番号を指定
retsu=4 #列番号を指定
xmax = system("cat " . loadname . " | awk \'NR==" . gyou ."{print $" . retsu . "}\'") # 特定の行列から値を拾う
set arrow 1 from xmax,y_min to xmax,y_max nohead lc rgb "forest-green" lw 1

f(x)=real(t1)
set view 0,0,,
spl loadname index 1 u (sqrt(($1)*($1)+($2)*($2)+($3)*($3))):4:5 w p pt 7 ps 0.3 palette notitle,\
