set terminal postscript; set output "../plots/plot10.ps"; set yrange [0:100]; set xlabel "year"; set ylabel "%"; set key left top;
plot "../data/results1.dat" u 1:(100*$9) title "L in HIV+" w li ,"../data/results1.dat" u 1:(100*$10) title "L in HIV-" w li,"../data/results3.dat" u 1:(100*$5) title "recent in HIV+" w li;
 #  "../data/pop.dat" u 1:2 poittype 7, "../data/TBnotes.dat" u 1:2 pointtye 5;\
# t,fast,smr+prev,smr+detections
# "../data/results2.dat" u 1:(100*$4) title "smr+notif" w li,#this was very noisy