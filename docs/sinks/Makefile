FNAME=SinksInTorch

all: ${FNAME}.pdf

%.pdf: %.tex Makefile
	pdflatex $<
	pdflatex $<

.PHONY: clean

clean:
	rm -f ${FNAME}.aux
	rm -f ${FNAME}.log
	rm -f ${FNAME}.out
	rm -f ${FNAME}.pdf
#rm ${FNAME}.bbl
#rm ${FNAME}.blg
#rm ${FNAME}Notes.bib
