#! /usr/bin/gnuplot
#filename = "file1.txt"

line_y_up = "zup.txt"

set xlabel "X-coord (M)"
set ylabel "Path (m)"
set grid

plot filename using 1:3 with lines title "z-path" lw 3 lt rgb "red",\
     line_y_up using 1:2 with lines title "up window" lw 5 lt rgb "blue",\
     line_y_up using 3:4 with lines title " down window" lw 5 lt rgb "blue"


