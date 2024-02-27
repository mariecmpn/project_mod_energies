# Exemple de script Gnuplot pour les solutions b du probleme inverse

set xlabel "t"
set ylabel "a(t)"
set title "Solutions a du problème inverse"
plot "solution_b.dat" with lines title "solution b_{app}", [0:1] 1+2*x with lines title "solution b_{ex}"