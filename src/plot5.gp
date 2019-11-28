# set terminal postscript; set output "../plots/plot05.ps"; set yrange [0:*]; set y2range [0:*]; set ytics nomirror; set xlabel "year"; set ylabel "% (interventions)"; set y2label "% (general ART coverage)"; set key left top;
# plot "../data/results5.dat" u 1:(100*$2) title "prevalence of HH" w li ax x1y1,"../data/results5.dat" u 1:(100*$3) title "prevalence of past HH ART" w li ax x1y1,"../data/results5.dat" u 1:(100*$4) title "prevalence of HH IPT" w li ax x1y1, "../data/results5.dat" u 1:(100*$7) title "prevalence ART in HIV+" w li ax x1y2,"../data/results5.dat" u 1:(100*$8) title "prevalence of ART in eligible" w li ax x1y2; #end of this line, start of plot
set terminal postscript;
set output "../plots/plot05.ps";
set size 1.0,0.5;
set multiplot; set origin 0.0,0.0; set size 0.5,1.0; set yrange [0:*]; set xrange [2000:*]; set xlabel "year"; set ylabel "% (interventions)"; plot "../data/results5.dat" u 1:(100*$2) title "prevalence of HH" w li ,"../data/results5.dat" u 1:(100*$3) title "prevalence of past HH IPT" w li ,"../data/results5.dat" u 1:(100*$4) title "prevalence of HH art" w li, "../data/results5.dat" u 1:(100*$5) title "prevalence of touched" w li ;
set origin 0.5,0.0; set size 0.5,1.0; set xrange [2000:*];
set ylabel "% (general ART coverage)"; set key left top;
plot "../data/results5.dat" u 1:(100*$7) title "ART in HIV+" w li ,"../data/results5.dat" u 1:(100*$8) title "ART in eligible" w li,"../data/results5.dat" u 1:(100*$15) title "HIV+ linked to care" w li ;
unset multiplot;
