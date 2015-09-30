set terminal epslatex color
set size ratio 1/1.618
set output 'irform.eps'
set xrange[3600:600]
set yrange[1:0]
set format xy '\small $%g$'
set xlabel 'Wavenumber (cm\textsuperscript{-1})'; set ylabel 'Relative Absorbance'
set key inside r b
set xtics 600, 400, 3600
plot 'form_ir.dat' u 1:2 w l lw 3 title '\footnotesize Simulated', 'form_ir.dat' u 3:(1-$4) w l lt 3 lc 3 lw 3 title '\footnotesize Experimental'

set size ratio 1/1.618
set output 'irtno3.eps'
set xrange[1200:0]
set yrange[1:0]
set format xy '$%g$'
set xlabel 'Wavenumber (cm\textsuperscript{-1})'; set ylabel 'Relative Absorbance'
set xtics 0,200,1200
plot 'no3_f_orca.txt' u 1:2 w l lw 3 lc 7 notitle
