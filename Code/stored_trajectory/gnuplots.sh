#!/bin/bash
for (( i = 0 ; i <= 250 ; i += 1 )) ; do
  gnuplot -e "filename='file$i.txt'" -p plot_z.gp 
done