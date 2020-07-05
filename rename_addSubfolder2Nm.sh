for d in fold*/sub*/; do file=`ls $d*.txt`; f=${d/\/*/}; echo $f; echo $file; outfile=${file/.txt/.$f.txt}; echo $outfile; mv $file $outfile; done
