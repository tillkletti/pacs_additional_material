TARGET=slides

AUXFILES=beamerthemePOLIMI.sty

PDFLATEX=pdflatex -synctex=1

$(TARGET).pdf: $(TARGET).tex $(AUXFILES)
	$(PDFLATEX) $<
	$(PDFLATEX) $<
	make clean


.PHONY: pack clean pdfclean

pack:
	tar czvf $(TARGET).tar.gz $(TARGET).tex $(AUXFILES) figures

clean:
	$(RM) $(TARGET).log $(TARGET).aux $(TARGET).out $(TARGET).snm $(TARGET).nav $(TARGET).toc *~ $(TARGET).syntex.gz

pdfclean: clean
	$(RM) $(TARGET).pdf
