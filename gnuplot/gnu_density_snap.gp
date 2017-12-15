reset
set size square
set ticslevel 0
set grid

set xr[0.1:10]
set yr[0.0001:10]
set logscale xy

set xl 'radius[kpc]' font 'Arial,15'
set yl 'density' font 'Arial,15'
filename='./data_0_0/density_0.dat'

set term png
set output './data_0_0/density_snap.png'

f(x)=(3/(4*pi))/(sqrt((1+(x**2)))**5)
pl filename u 1:2 with boxes notitle,f(x) w l title 'Analytical'
