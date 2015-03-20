.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:bfc.pdf supp.pdf

bfc.pdf:bfc.tex bfc.bib
		pdflatex bfc; bibtex bfc; pdflatex bfc; pdflatex bfc;

supp.pdf:supp.tex bfc.bib
		pdflatex supp; bibtex supp; pdflatex supp; pdflatex supp;

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~
