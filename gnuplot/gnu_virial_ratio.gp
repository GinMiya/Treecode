reset
set xr[0:200]
set yr[0:1.5]
set ticslevel 0
set view equal xy
set size square
set grid
set term png
dirname = './data_0_0'
filename = dirname.'/virial_ratio.png'
loadname = dirname.'/energy_data.dat'
set output filename
set xl 'time' font 'Arial,20'
set yl 'Virial Ratio' font 'Arial,20'
f(x)=1.0
pl loadname u 1:4 smooth csplines w l title 'Virial'\
,f(x) w l title '1.0'\
,loadname u 1:2 smooth csplines w l title 'kinetic'\
,loadname u 1:(-$3) smooth csplines w l title '|potential|'\
,loadname u 1:6 smooth csplines w l title 'momentum'
