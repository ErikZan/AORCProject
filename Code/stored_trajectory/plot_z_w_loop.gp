#! /usr/bin/gnuplot
#filename = "file1.txt"

line_y_up = "zup.txt"

set xlabel "X-coord (M)"
set ylabel "Path (m)"
set grid

do for [i=1:250] {
    set output sprintf('ftg/%d.png', i)
    filename=sprintf('file%d.txt', i)
    plot filename using 1:3 with lines title "y-path" lw 3 lt rgb "red",\
    line_y_up using 1:2 with lines title "left window" lw 5 lt rgb "blue",\
    line_y_up using 3:4 with lines title " right window" lw 5 lt rgb "blue"

}




