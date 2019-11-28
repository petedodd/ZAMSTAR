set terminal postscript; set output "../plots/plot04b.ps"; set yrange [0:*]; set xlabel "year"; set ylabel "per 100,000 per year"; set key right top;
plot  "../data/results4.dat" u 1:5 title "TB incidence" w li , "../shared/GDsainc.dat" u 1:2:3:4 ti "SA incidence" w yerr ;
