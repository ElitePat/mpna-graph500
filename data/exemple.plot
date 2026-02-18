set title "Trois courbes avec gnuplot"
set xlabel "X"
set ylabel "Y"
set grid

plot "data.txt" using 1:2 with lines title "Courbe 1", \
     "data.txt" using 1:3 with lines title "Courbe 2", \
     "data.txt" using 1:4 with lines title "Courbe 3"
