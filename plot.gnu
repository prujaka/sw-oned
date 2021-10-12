set terminal postscript eps enhanced defaultplex color size 10cm,15cm \
font 'Courier,16'
set output 'hu.eps'
set style data lines
set linestyle 1 lw 3.0 lc rgb 'red'
set linestyle 2 lw 3.0 lc rgb 'blue'
filename = 'res.dat'

set key top left
#set yrange[0.9:1.9]

set xlabel 'x, m'

set multiplot layout 2,1

set ylabel 'h, m'
plot filename u 1:2 ls 1 title ''

set ylabel 'u, m'
plot filename u 1:3 ls 1 title ''

unset multiplot
