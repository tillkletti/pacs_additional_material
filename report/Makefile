
MAIN = main
# Nom du fichier principal .tex

TITLE = report
# Nom souhaité du fichier final .pdf


FOLDER = report

default:
	make draft

$(TITLE).pdf: $(MAIN).pdf
	rm -f $(TITLE).pdf
	mv $(MAIN).pdf $(TITLE).pdf
	make distclean
	

$(MAIN).pdf: $(MAIN).tex
	pdflatex $(MAIN)
	#pdflatex $(MAIN)
	#makeglossaries $(MAIN)
	bibtex $(MAIN)
	pdflatex $(MAIN)
	#pdflatex $(MAIN)

draft: $(MAIN).tex
	pdflatex $(MAIN)
	make clean

proper:
	make $(TITLE).pdf

clean: 
	rm -f $(MAIN).log
	rm -f $(MAIN).aux
	rm -f $(MAIN).alg
	rm -f $(MAIN).acr
	rm -f $(MAIN).acn
	rm -f $(MAIN).bbl
	rm -f $(MAIN).bcf
	rm -f $(MAIN).blg
	rm -f $(MAIN).lof
	rm -f $(MAIN).out
	#rm -f $(MAIN).pdf
	rm -f $(MAIN).toc
	rm -f $(MAIN).run.xml
	rm -f $(MAIN).synctex.gz
	rm -f $(MAIN).glo
	rm -f $(MAIN).gls
	rm -f $(MAIN).glg
	rm -f $(MAIN).ist
	rm -f $(MAIN)-blx.bib
	rm -f *.aux
	rm -f *.log

distclean:
	make clean
	rm -f $(MAIN).pdf

