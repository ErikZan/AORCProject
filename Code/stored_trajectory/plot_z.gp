#! /usr/bin/gnuplot
#filename = "file1.txt"

set xlabel "X-coord (M)"
set ylabel "Path (m)"
set grid


plot filename using 2:4 with lines title "z-path"



