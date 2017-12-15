reset
set xr[-20:20]
set yr[-20:20]
set ticslevel 0
set view equal xy
set size square
set grid
set term gif animate delay 5
dirname = "./data_0_0"
filename = dirname."/position2d.gif"
set output filename
unset title

n0=0
nm=5000
dn=25

load 'gnu_position2d_gif.plt'
