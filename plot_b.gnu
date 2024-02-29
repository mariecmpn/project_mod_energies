# Exemple de script Gnuplot pour les solutions b du probleme inverse

set xlabel "t"
set ylabel "a(t)"
set title "Solutions b du problème inverse"
#set title "Graphique d'erreur pour b"
plot "solution_b.dat" with lines title "solution b_{app}", [0:1] 1+2*x with lines title "solution b_{ex}"
#plot "erreur_b.dat" with lines title "|b_{app} - b_{ex}|"