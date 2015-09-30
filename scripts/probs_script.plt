set terminal eps color enhanced#epslatex color enhanced
set output 'results.eps'
set xrange[0:30]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
plot 'probabilities.dat' u 1:2 w l title '$R_{1 \leftarrow 1}$', 'probabilities.dat' u 1:3 w l title '$R_{2 \leftarrow 1}$', 'probabilities.dat' u 1:4 w l title '$T_{1 \leftarrow 1}$', 'probabilities.dat' u 1:5 w l title '$T_{2 \leftarrow 1}$'
