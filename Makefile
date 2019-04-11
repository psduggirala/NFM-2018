
MAIN=emsoft-2019

all $(MAIN).pdf : $(MAIN).tex 
	pdflatex $<
	-bibtex $(basename $<)
	pdflatex $<
	pdflatex $<

STEM=$(notdir $(PWD))

snapshot : $(MAIN).pdf
	@TODAY="`date +%Y%m%d`"; \
	N="1"; \
	while [ -f $$TODAY-$(STEM)-$$N.pdf ]; do \
	  N="`expr $$N + 1`"; done; \
	echo $(MAIN).pdf '->' $$TODAY-$(STEM)-$$N.pdf; \
	cp $(MAIN).pdf $$TODAY-$(STEM)-$$N.pdf

%.pdf : %.gpl
	gnuplot $(<)

plots : $(subst gpl,pdf,$(wildcard *.gpl))

clean :
	rm -f $(addprefix $(MAIN)., aux bbl blg fmt log out pdf)

.PHONY : plots snapshot clean

