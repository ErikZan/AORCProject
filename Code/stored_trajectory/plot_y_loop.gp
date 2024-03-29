#! /usr/bin/gnuplot

set xlabel "X-coord (M)"
set ylabel "Path (m)"
set grid
set terminal pngcairo
set pm3d map

do for [i=1:250] {
    set output sprintf('ftg/%d.png', i)
    filename=sprintf('file%d.txt', i)
    plot filename using 1:2 with lines title "y-path"
}