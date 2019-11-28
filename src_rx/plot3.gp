set terminal postscript; set output "../plots/plot03.ps"; set yrange [0:*] ;  set y2range [0:*];set ytics nomirror; set y2tics nomirror; set xlabel "year"; set ylabel "% per year"; set y2label "%"; set key left top;
plot "../data/results3.dat" u 1:(100*$2) title "ARI>12" w li ax x1y1,"../data/results3.dat" u 1:(100*$3) title "ARI<12" w li ax x1y1,"../data/results3.dat" u 1:(100*$4) title "TB prev HH" w poi ax x1y2;
 #  "../data/pop.dat" u 1:2 poittype 7, "../data/TBnotes.dat" u 1:2 pointtye 5;
#t,FOI,ARI