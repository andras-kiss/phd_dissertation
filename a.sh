# gfortran ULL_deconvolution.F90 -o ULL_a.out
# ./ULL_a.out
# gnuplot ULL_plot

#gnuplot lineplot
#ls *.eps > eps_files
#for i in $(cat eps_files)
#do
#   ps2pdf $i
#done
#gnuplot plot_transient

# gnuplot plot_sb
pdflatex phd.tex
bibtex phd
pdflatex phd.tex
pdflatex phd.tex

rm *.{aux,bbl,blg,lof,log,lot,out,toc}
