# makefile for thesis.pdf
filename=thesis

pdf: thesis.pdf
	bibtex ${filename}; pdflatex -interaction=nonstopmode ${filename}.tex
	bibtex ${filename}; pdflatex -interaction=nonstopmode ${filename}.tex
view:
	evince thesis.pdf &

cleanall:
	 rm -f ${filename}.log ${filename}.aux ${filename}.out ${filename}.lot ${filename}.bbl ${filename}.blg ${filename}.bfr ${filename}.glo ${filename}.toc
	 rm -f ${filename}.ps ${filename}.pdf 

clean:
	 rm -f ${filename}.log ${filename}.aux ${filename}.out ${filename}.lot ${filename}.bbl ${filename}.blg ${filename}.brf ${filename}.glo ${filename}.toc
