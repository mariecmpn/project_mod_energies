# Exemple de script Gnuplot pour afficher la solution du pb inverse

set xlabel "x"
set ylabel "t"
set zlabel "u(x,t)"
splot "solution_u.dat" u 1:2:3 title "solution approchée", [0:1] [0:1] (exp(x)+x**2)*exp(-y) title "solution exacte"

