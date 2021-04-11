# Gnuplot script file
# Automatically generated by eggplot Ver. 0.1.0
set datafile separator ','
set terminal pngcairo dashed enhanced
set output 'eggp-export.png'
set grid lc rgb '#cccccc' lw 1 lt 3
set style line 1 lt 1 lw 1 pt 1 ps 1 lc rgb '#f00032'
set style line 2 lt 1 lw 1 pt 2 ps 1 lc rgb '#227500'
set style increment userstyle
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title ""
set xlabel "{/Symbol D}P"
set ylabel "v_0"
plot 'eggp.dat' index 0 title '1' with linespoints, 'eggp.dat' index 1 title '2' with linespoints, 
