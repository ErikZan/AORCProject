line_y_up = "yup.txt"
line_z_up = "zup.txt"

set xlabel "X-coord (M)"
set ylabel "Path y (m)"
set zlabel "Path z (m)"

set grid
set parametric;

splot filename using 1:2:3 with lines title "z-path" lw 3 lt rgb "red",\
     line_y_up using 1:2:(1.5) with lines title "left window" lw 5 lt rgb "blue",\
     line_y_up using 3:4:(2.5) with lines title " right window" lw 5 lt rgb "blue",\
     line_z_up using 1:2:(0.0) with lines title "down window" lw 5 lt rgb "blue",\
     line_z_up using 3:4:(0.0) with lines title " up window" lw 5 lt rgb "blue",\