# Exemple de script Gnuplot pour les solutions a du probleme inverse

set xlabel "t"
set ylabel "a(t)"
set title "Solutions a du problème inverse"
plot "solution_a.dat" with lines title "solution a_{app}", [0:1] 1+x with lines title "solution a_{ex}"