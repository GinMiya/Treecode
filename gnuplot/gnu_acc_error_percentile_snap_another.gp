reset
set ticslevel 0
set size square
set grid
#set logscale y
set yr[-10**(-5):10**(5)]
set xr[0:100]
set key left
set term png
dirname = '../data_0_0'
#dirname = '../cold_collapse'
filename = dirname.'/acc_error_percentile.png'
loadname = dirname.'/acc_error.dat'
set output filename
set format y "10^{%L}"
set xl 'Force error magnitude' font "Arial,15"
set yl 'Percentile' font 'Arial,15'
plot loadname u 1:2 w l notitle
