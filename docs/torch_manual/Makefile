all: main.pdf

main.pdf: main.tex Makefile fig/pe_heating.pdf fig/gas_cool.pdf fig/recombination_case_b.pdf fig/dust_opacity.pdf fig/dust_heat_cool_Tgas10.pdf fig/dust_heat_cool_Tgas10000.pdf fig/gas_atomic_cool.png fig/strong_scaling_cart.png fig/weak_scaling_cart.png
	pdflatex main.tex
	pdflatex main.tex

fig/pe_heating.pdf: pe_heating.py
	mkdir -pv fig
	./pe_heating.py

fig/gas_cool.pdf: gas_cool.py cool.dat
	mkdir -pv fig
	./gas_cool.py

fig/recombination_case_b.pdf: recombination_case_b.py
	mkdir -pv fig
	./recombination_case_b.py

fig/dust_opacity.pdf: dust_opacity.py
	mkdir -pv fig
	./dust_opacity.py

fig/dust_heat_cool_Tgas10.pdf fig/dust_heat_cool_Tgas10000.pdf: dust_heat_cool.py
	mkdir -pv fig
	./dust_heat_cool.py

.PHONY: clean

# this should probably clear out the figures, too
clean:
	rm -f main.aux
	rm -f main.log
	rm -f main.out
	rm -f main.pdf
	rm -f main.toc
