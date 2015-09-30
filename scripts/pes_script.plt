set terminal epslatex color# enhanced

set output 'pes_sc.eps'
set xrange[-4:4]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
plot 'sac.dat' u 1:2 w l lc 1 lw 3 title '$H_{11}$', 'sac.dat' u 1:3 w l lt 4 lc 4 lw 3 title '$H_{12}=H_{21}$', 'sac.dat' u 1:5 w l lt 2 lc 3 lw 3 title '$H_{22}$'

set output 'a_pes_sc.eps'
set xrange[-4:4]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
plot 'sac.dat' u 1:6 w l lw 3 title '$E_{1}$', 'sac.dat' u 1:7 w l lc 3 lw 3 title '$E_{2}$

set output 'dpes_sc.eps'
set xrange[-4:4]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
set key spacing 1.75
plot 'dsac.dat' u 1:2 w l lc 1 lw 3 title '$\dpar{H_{11}}{R}$', 'dsac.dat' u 1:3 w l lt 4 lc 4 lw 3 title '$\dpar{H_{12}}{R} = \dpar{H_{21}}{R}$', 'dsac.dat' u 1:5 w l lt 2 lc 3 lw 3 title '$\dpar{H_{22}}{R}$'

set output 'pes_dc.eps'
set xrange[-8:8]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
plot 'dac.dat' u 1:2 w l lc 1 lw 3 title '$H_{11}$', 'dac.dat' u 1:3 w l lt 4 lc 4 lw 3 title '$H_{12}=H_{21}$', 'dac.dat' u 1:5 w l lt 2 lc 3 lw 3 title '$H_{22}$'

set output 'dpes_dc.eps'
set xrange[-8:8]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
set key spacing 1.75
plot 'ddac.dat' u 1:2 w l lc 1 lw 3 title '$\dpar{H_{11}}{R}$', 'ddac.dat' u 1:3 w l lt 4 lc 4 lw 3 title '$\dpar{H_{12}}{R} = \dpar{H_{21}}{R}$', 'ddac.dat' u 1:5 w l lt 2 lc 3 lw 3 title '$\dpar{H_{22}}{R}$'

set output 'a_pes_dc.eps'
set xrange[-8:8]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
plot 'dac.dat' u 1:6 w l lw 3 title '$E_{1}$', 'dac.dat' u 1:7 w l lc 3 lw 3 title '$E_{2}$

set output 'pes_ec.eps'
set xrange[-10:10]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
plot 'ec.dat' u 1:2 w l lc 1 lw 3 title '$H_{11} \times 50$', 'ec.dat' u 1:3 w l lt 4 lc 4 lw 3 title '$H_{12}=H_{21}$', 'ec.dat' u 1:5 w l lt 2 lc 3 lw 3 title '$H_{22} \times 50$'

set output 'dpes_ec.eps'
set xrange[-10:10]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
set key spacing 1.75
plot 'dec.dat' u 1:2 w l lc 1 lw 3 title '$\dpar{H_{11}}{R} = \dpar{H_{22}}{R}$', 'dec.dat' u 1:3 w l lt 4 lc 4 lw 3 title '$\dpar{H_{12}}{R} = \dpar{H_{21}}{R}$'

set output 'a_pes_ec.eps'
set xrange[-10:10]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
plot 'ec.dat' u 1:6 w l lw 3 title '$E_{1}$', 'ec.dat' u 1:7 w l lc 3 lw 3 title '$E_{2}$

set output 'da_pes_ec.eps'
set xrange[-10:10]
set format xy '$%g$'
set grid
set xlabel 'Position ($\mathbf{R}$)'; set ylabel 'Energy (A.U.)'
set key spacing 1.75
plot 'dec.dat' u 1:6 w l lw 3 title '$\dpar{E_{1}}{R}$', 'dec.dat' u 1:7 w l lc 3 lw 3 title '$\dpar{E_{2}}{R}$
