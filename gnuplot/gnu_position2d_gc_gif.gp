reset
set xr[-50:50]
set yr[-50:50]
set ticslevel 0
set view equal xy
set size square
set grid
set term gif animate delay 5
dirname = "./data_0_0"
filename = dirname."/position2d_gc.gif"
set output filename
unset title

n0=0
nm=5000
dn=25

load 'gnu_position2d_gc_gif.plt'
