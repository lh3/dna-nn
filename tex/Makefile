.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

.pdf.eps:
		pdftops -eps $< $@

all:dna-nn.pdf

dna-nn.pdf:dna-nn.tex dna-nn.bib
		pdflatex dna-nn; bibtex dna-nn; pdflatex dna-nn; pdflatex dna-nn;

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ dna-nn.pdf
