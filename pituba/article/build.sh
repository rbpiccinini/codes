#!/bin/bash

pdflatex article_spe.tex
bibtex article_spe
pdflatex article_spe.tex
pdflatex article_spe.tex

evince article_spe.pdf &
exit 1

