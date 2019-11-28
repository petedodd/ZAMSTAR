set terminal postscript; set output "../plots/plot02.ps"; set yrange [0:100]; set xlabel "year"; set ylabel "%"; set key left top;
plot "../data/results2.dat" u 1:(100*$2) title "recent" w li ,"../data/results2.dat" u 1:(100*$3) title "smr+prev" w li ,"../data/results2.dat" u 1:(100*$5) title "FOI via HH" w li, "../data/results2.dat" u 1:(100*$6) w p pt 2 title "HIV in HH1549", "../data/results2.dat" u 1:(100*$8) w p pt 6 title "Recurrent in Rx" ;
 #  "../data/pop.dat" u 1:2 poittype 7, "../data/TBnotes.dat" u 1:2 pointtye 5;\
# t,fast,smr+prev,smr+detections
# "../data/results2.dat" u 1:(100*$4) title "smr+notif" w li,#this was very noisy