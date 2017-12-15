reset
set xr[0:200]
set yr[0.00001:0.01]
set ticslevel 0
set view equal xy
set size square
set logscale y
set grid
set term png
dirname = './data_0_0'
filename = dirname.'/energy_error.png'
loadname = dirname.'/energy_data.dat'
set output filename
set xl 'time' font 'Arial,20'
set yl 'energy error' font 'Arial,20'
pl loadname every ::2 u 1:5 smooth acsplines w l title 'Energy error'
