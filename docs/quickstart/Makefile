all: quickstart-web.pdf
#quickstart-cca.pdf

quickstart-cca.pdf: quickstart-cca.tex fig/torch-features.pdf fig/quickstart-particles.pdf fig/quickstart-dens-slice.pdf
	pdflatex quickstart-cca.tex
	pdflatex quickstart-cca.tex

quickstart-web.pdf: quickstart-web.tex fig/torch-features.pdf fig/quickstart-particles.pdf fig/quickstart-dens-slice.pdf
	pdflatex quickstart-web.tex
	pdflatex quickstart-web.tex

.PHONY: clean

clean:
	rm -f quickstart-web.aux
	rm -f quickstart-web.log
	rm -f quickstart-web.out
	rm -f quickstart-web.pdf
	rm -f quickstart-web.toc
