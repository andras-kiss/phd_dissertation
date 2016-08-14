gfortran deconvolution.F90
./a.out
#gnuplot plot
#./autoplot.sh

gnuplot lineplot

ls *.eps > eps_files
for i in $(cat eps_files)
do
   ps2pdf $i
done
