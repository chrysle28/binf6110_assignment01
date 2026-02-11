set terminal png size 2000,2000
set output "dotplot.png"

set xlabel "Reference genome (bp)"
set ylabel "Assembly contigs (bp)"

set grid
set key off

plot "dotplot.coords" using 1:2:($3-$1):($4-$2) with vectors nohead lw 2 lc rgb "blue"
