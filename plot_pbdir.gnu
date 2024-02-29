# Exemple de script Gnuplot pour afficher la solution du pb direct

set xlabel "x"
set ylabel "t"
set zlabel "u(x,t)"
set title "Solutions du problème direct"
#set title "Graphique d'erreur pour u (problème direct)"
splot "solution_u_pbdir.dat" u 1:2:3 title "solution approchée", [0:1] [0:1] (exp(x)+x**2)*exp(-y) title "solution exacte"
#splot [0:1] [0:1] (exp(x)+x**2)*exp(-y) title "solution exacte"
#splot "erreur_u_pbdir.dat" u 1:2:3 title "|u_{app}-u_{ex}|"