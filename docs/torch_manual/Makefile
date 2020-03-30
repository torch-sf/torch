all: main.pdf quickstart-web.pdf
#quickstart-cca.pdf

quickstart-cca.pdf: quickstart-cca.tex torch-features.pdf quickstart-particles.pdf quickstart-dens-slice.pdf
	pdflatex quickstart-cca.tex
	pdflatex quickstart-cca.tex

quickstart-web.pdf: quickstart-web.tex torch-features.pdf quickstart-particles.pdf quickstart-dens-slice.pdf
	pdflatex quickstart-web.tex
	pdflatex quickstart-web.tex

main.pdf: main.tex fig/pe_heating.pdf fig/gas_cool.pdf fig/recombination_case_b.pdf fig/dust_opacity.pdf fig/dust_heat_cool_Tgas10.pdf fig/dust_heat_cool_Tgas10000.pdf fig/pe_heating_expressions.pdf fig/gas_atomic_cool.png
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

fig/pe_heating_expressions.pdf: pe_heating_expressions.py
	mkdir -pv fig
	./pe_heating_expressions.py

# TOTAL hack -AT 2019 June 21
# IDK better way to do this without having to git commit binary files...
fig/gas_atomic_cool.png:
	mkdir -pv fig
	wget -O fig/gas_atomic_cool.png http://user.astro.columbia.edu/~atran/share/gas_atomic_cool.png
