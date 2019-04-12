set terminal pdfcairo enhanced dashed font 'Verdana,18' fontscale 0.45 dl 0.4

set size ratio 0.5
set lmargin at screen 0.12
set rmargin at screen 0.99
set bmargin at screen 0.15
set tmargin at screen 0.9

set autoscale
set yrange [0.05:65]
set xrange [0:72]
unset label
set xtic auto
set ytic auto
set xlabel "Passive Switch Time Radius (r)" offset 0,0.5
set ylabel "Runtime (seconds)" offset 2.5,0

set style line 1 lc rgb "#00bd60" lt 1 lw 3.5 pt 12 ps 0.5
set style line 2 lc rgb "#ff4040" lt 2 lw 3.0 pt 5 ps 0.2 dt 3
set style line 3 lc rgb "#8B008B" lt 6 lw 0.6 pt 3 ps 0.4

set pointintervalbox 0.5
set style line 4 lc rgb "#ef8d00" lt 2 lw 2.0 pt 6 ps 0.4 dt 2 pi -1

set key on bottom right
set logscale y

set output "hylaa_unagg.pdf"
set title "Hylaa without Aggregation" font 'Verdana,22' offset 0,-1

set style line 5 lc rgb "#808080" lt 2 lw 2.0 pt 13 ps 0.8 dt 2
set arrow from 0, 50 to 75, 50 nohead ls 5
set label 21 "Timeout (One Minute)" at 10,50 font 'Verdana,14' tc rgb "#808080" offset 0, -0.4

set grid ytics lc '#606060' lw 0.25 lt 1 dt 3
set grid xtics lc '#606060' lw 0.25 lt 1 dt 3

set key Left invert font 'Verdana,14' box

datafile = 'data_hylaa_unagg.dat'
plot for [i=0:*] datafile index i using 1:2\
with linespoints title columnheader(1) ls (i+1)
