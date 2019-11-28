set terminal postscript; set output "../plots/plot06.ps"; set yrange [0:*]; set xlabel "year"; set ylabel "%"; set key left top;
plot "../data/results5.dat" u 1:(100*$5) title "prevalence of ECF" w li,"../data/results5.dat" u 1:(100*$6) title "ECF detection proportion" w li;
 #  "data/pop.dat" u 1:2 poittype 7, "data/TBnotes.dat" u 1:2 pointtye 5;\
# t,fast,smr+prev,smr+detections
# "data/results5.dat" u 1:(100*$4) title "smr+notif" w li,#this was very noisy