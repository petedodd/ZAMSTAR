set terminal postscript;
 set output "../plots/plot15.ps"; 
 set yrange [0:*]; 
 set y2range [0:*];
 set ytics nomirror; 
 set y2tics nomirror; 
 set xlabel "year"; set ylabel "per 100,000"; set y2label "per 100,000 per year"; set key right top;
plot "results6.dat" u 1:2 title "Prev on standard Rx (1st line)" w li axis x1y1 , "results6.dat" u 1:3 title "Prev on REMox Rx (1st line)" w li axis x1y1, "results6.dat" u 1:4 title "Prev on standard Rx (retreatment)" w li axis x1y1, "results6.dat" u 1:5 title "Prev on REMox Rx (retreatment)" w li axis x1y1, "results6.dat" u 1:6 title "Number starting Standard" w li axis x1y2, "results6.dat" u 1:7 title "Number starting REMox" w li axis x1y2 ;
