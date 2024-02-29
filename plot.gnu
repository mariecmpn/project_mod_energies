# Exemple de script Gnuplot pour afficher la solution du pb inverse

set xlabel "x"
set ylabel "t"
set zlabel "u(x,t)"
set title "Solutions u du problème inverse"
#set title "Graphique d'erreur pour u"
splot "solution_u.dat" u 1:2:3 title "solution approchée", [0:1] [0:1] (exp(x)+x**2)*exp(-y) title "solution exacte"
#splot "erreur_u.dat" u 1:2:3 title "|u_{ex} - u_{app} |"

