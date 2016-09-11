#!/usr/bin/gnuplot

set logscale y
set title "Residuals"
set ylabel 'Residual'
set xlabel 'Iteration'
plot " < tail -n 500 log2 | grep 'Liquid Mass in system' | cut -d' ' -f6" title 'liquid mass' with lines 
pause -1
