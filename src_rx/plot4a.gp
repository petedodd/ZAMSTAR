set terminal postscript; set output "../plots/plot04a.ps"; set yrange [0:*];  set xlabel "year"; set ylabel "per 100,000"; set key right top;
plot "../data/results4.dat" u 1:2 title "TB prevalence (18+)" w li  ,  "../shared/GDsaprev.dat" u 1:2:3:4 ti "SA prevalence" w yerr ;
