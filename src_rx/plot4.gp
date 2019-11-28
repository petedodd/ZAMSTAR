set terminal postscript; set output "../plots/plot04.ps"; set yrange [0:*]; set y2range [0:*];set ytics nomirror; set y2tics nomirror; set xlabel "year"; set ylabel "per 100,000"; set y2label "per 100,000 per year"; set key right top;
plot "../data/results4.dat" u 1:2 title "TB prevalence (18+)" w li axis x1y1 , "../data/results4.dat" u 1:5 title "TB incidence" w li axis x1y2, "../data/results4.dat" u 1:4 title "TB Rx (18+)" w li axis x1y1, "../data/results4.dat" u 1:3 title "TB Rx incidence" w li axis x1y2 ;
#, "../data/GDsainc.dat" u 1:2:3:4 ti "SA incidence" w yerrorlines axis x1y2, "../data/GDsaprev.dat" u 1:2:3:4 ti "SA prevalence" w yerrorlines axis x1y1;
#
 #  "../data/pop.dat" u 1:2 poittype 7, "../data/TBnotes.dat" u 1:2 pointtye 5;
#t,TB prev, notifs, Rx,Inc