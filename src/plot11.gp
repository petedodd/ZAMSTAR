set terminal postscript;
set output "../plots/plot11.ps";
set xlabel "year"; set ylabel "Coverage (%)";set key left top;
plot "../data/results5.dat" u 1:(100*$16) title "prevalence IPT" w li ,"../data/results5.dat" u 1:(100*$17) title "prevalence IPT, HIV+" w li,"../data/results5.dat" u 1:(100*$20) title "ever IPT among HIV+" w li;;
#"../data/results5.dat" u 1:(100*$18) title "prevalence IPT, ART" w li,