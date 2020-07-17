#! /usr/bin/gnuplot
#filename = "file1.txt"

set xlabel "X-coord (M)"
set ylabel "Path y (m)"
set zlabel "Path z (m)"

set grid

#plot filename using 2:4 with lines title "y-path"
#plot filename using 2:4 with lines title "z-path"
#splot filename using 2:3:4 with lines title "trajectory"
splot filename using 1:2:3 with lines title "trajectory"
