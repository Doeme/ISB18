
set terminal pdfcairo
set output "2_a.pdf"
set xrange [0:1]
set xlabel "V_{bio}"
set ylabel "V_{efni}"
plot (1-x)
set output "2_b.pdf"
plot (x)
