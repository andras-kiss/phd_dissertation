#!/bin/bash
ls *.txt > files
for i in $(cat files)
do
   gnuplot -e "filename='$i'; outputname='${i%.*}.eps'" autoplot
done

ls *diff.txt > files
for i in $(cat files)
do
   gnuplot -e "filename='$i'; outputname='${i%.*}.eps'" autoplot_diff
done

ls *.eps > files
for i in $(cat files)
do
   ps2pdf $i
done
