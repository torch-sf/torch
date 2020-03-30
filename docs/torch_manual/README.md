README

To download this repository and build the quickstart guide AND Aaron's "Torch
Manual" (i.e., personal disorganized notes on Torch), do:

    git clone git@bitbucket.org:aarontran/torch_manual.git
    cd torch_manual
    make

To contribute: feel free to push directly to master, but I recommend to `git
pull` and make sure there aren't conflicting commits, before you push.

This private repository contains (as of Feb 26):

__Torch quickstart guide.__  Most recent version is "-web".  This is the
quickstart guide posted at:
[https://torch-sf.bitbucket.io/quickstart.pdf](https://torch-sf.bitbucket.io/quickstart.pdf).
Relevant files are:

	quickstart-cca.tex
	quickstart-web.tex
	quickstart-dens-slice.pdf
	quickstart-particles.pdf
	torch-features.pdf

__Aaron's personal notes on Torch ("A Torch Manual").__  Some overlap with
quickstart (basically, a subset of these notes evolved into the quickstart).
A copy is posted at:
[http://user.astro.columbia.edu/~atran/torch/manual.pdf](http://user.astro.columbia.edu/~atran/torch/manual.pdf).
Username/password: `torch`/`torch`.
Relevant files are:

	main.tex
	post_atran.sh

Scripts to make figures in aforementioned notes:

	cool.dat
	dust_heat_cool.py
	dust_opacity.py
	gas_cool.py
	pe_heating.py
	pe_heating_expressions.py
	recombination_case_b.py
