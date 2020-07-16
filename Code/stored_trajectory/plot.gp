#! /usr/bin/gnuplot
# file = "CutBoundary_Sensitivity.txt"

#file = "esperimento_gnuplot.txt"

#file=system("echo $name")
set xlabel "X-coord (M)"
set ylabel "Path (m)"

plot filename using 1:2 with lines
