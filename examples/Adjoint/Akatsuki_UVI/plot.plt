reset

set linetype 1 lc rgb '#1f77b4'
set linetype 2 lc rgb '#ff7f0e'
set linetype 3 lc rgb '#2ca02c'
set linetype 4 lc rgb '#d62728'
set linetype 5 lc rgb '#9467bd'
set linetype 6 lc rgb '#8c564b'
set linetype 7 lc rgb '#e377c2'
set linetype 8 lc rgb '#7f7f7f'
set linetype 9 lc rgb '#bcbd22'
set linetype 10 lc rgb '#17becf'
set linetype cycle 10

set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')

set term pngcairo size 1000, 1200 font "Open Sans, 20"

# set output "vertical_sensitivity_cloud.png"
# set xlabel "Sensitivity"
# set ylabel "Altitude (km)"
# set grid
# set ytics 5
# set mytics 5
# set key left bottom
# set logscale x
# set xrange[1e-9:1e15]
# set yrange[50:80]
# plot "vertical_sensitivity_cloud.dat" u 3:($1*1e-3) w l lw 3 title "{/Symbol s}_g (mode 1)", "vertical_sensitivity_cloud.dat" u 5:($1*1e-3) w l lw 3 title "{/Symbol s}_g (mode 2)", "vertical_sensitivity_cloud.dat" u 7:($1*1e-3) w l lw 3 title "{/Symbol s}_g (mode 2p)", "vertical_sensitivity_cloud.dat" u 9:($1*1e-3) w l lw 3 title "{/Symbol s}_g (mode 3)"

# set output "vertical_sensitivity_cloud.png"
# set xlabel "Sensitivity"
# set ylabel "Altitude (km)"
# set grid
# set ytics 5
# set mytics 5
# # set logscale x
# set key bottom
# # plot "vertical_sensitivity_cloud.dat" u 2:($1*1e-3) w l lw 3 title "r_g (mode 1)", "vertical_sensitivity_cloud.dat" u 4:($1*1e-3) w l lw 3 title "r_g (mode 2)", "vertical_sensitivity_cloud.dat" u 6:($1*1e-3) w l lw 3 title "r_g (mode 2p)", "vertical_sensitivity_cloud.dat" u 8:($1*1e-3) w l lw 3 title "r_g (mode 3)",\
# plot "vertical_sensitivity_cloud.dat" u 3

set output "vertical_sensitivity.png"
set xlabel "Relative Jacobian q_{SO_2}(dJ/dq_{SO_2})"
set ylabel "Altitude (km)"
set grid
set ytics 5
set mytics 5
set key bottom
plot "vertical_sensitivity.dat" u 7:($1*1e-3) w l lw 3 title "SO_2"


# set term pngcairo size 1100, 1200 font "Open Sans, 20"
# set output "vertical_phase_function.png"
# set xlabel "Scattering angle (degree)"
# set ylabel "Altitude (km)"
# set cblabel "dJ/dP({/Symbol q})"
# set grid
# set size ratio 1.5
# set logscale cb
# set ytics 5
# set xrange[0:180]
# set xtics 30
# set yrange[50:80]
# set cbrange[1e-20:1e-12]
# # set key bottom
# plot "vertical_phase_function.dat" u ($1*180/pi):($2*1e-3):3 w image notitle