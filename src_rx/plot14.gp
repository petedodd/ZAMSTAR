set terminal postscript; 
set output "../plots/plot14.ps"; 
set yrange [0:1]; 
set xlabel "year"; 
set ylabel "Proportion starting treatment that have been previously treated"; 
set key left top; 
set xrange [1980:2050];
plot "../data/results5.dat" u 1:($22) title "proportion" w li 
#
#, "../shared/popsg.dat" u 1:(100*($6-1)) w p pt 4 ps 2 lw 3 ti "ZMpop";
 #  "../data/pop.dat" u 1:2 poittype 7, "../data/TBnotes.dat" u 1:2 pointtye ;"../data/results1.dat" u 1:(100*$7) title "pop<12" w li,
 # t,deltaN,L,hiv,midage,hivinT
# "../data/results1.dat" u 1:(100*$11/$3) title "rL" w li,
