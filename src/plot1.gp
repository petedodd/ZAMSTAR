set terminal postscript; set output "../plots/plot01.ps"; set yrange [0:100]; set xlabel "year"; set ylabel "%"; set key left top; set xrange [1980:2050];
plot "../data/results1.dat" u 1:(100*$2) title "growth" w li ,"../data/results1.dat" u 1:(100*$3) title "L" w li ,"../data/results1.dat" u 1:(100*$4) title "HIV15-49" w li,"../data/results1.dat" u 1:(100*$5) title "pop15-49" w li, "../data/results1.dat" u 1:(100*$6) w p pt 2 ti "HIV TBRx", "../data/results2.dat" u 1:(100*$7) ti "HIV incTB" w li, "../shared/popsg.dat" u 1:(100*($5-1)) w p pt 3 ps 2 lw 3 ti "SApop";
#
#, "../shared/popsg.dat" u 1:(100*($6-1)) w p pt 4 ps 2 lw 3 ti "ZMpop";
 #  "../data/pop.dat" u 1:2 poittype 7, "../data/TBnotes.dat" u 1:2 pointtye ;"../data/results1.dat" u 1:(100*$7) title "pop<12" w li,
 # t,deltaN,L,hiv,midage,hivinT
# "../data/results1.dat" u 1:(100*$11/$3) title "rL" w li,