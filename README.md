README
======

Torch is an open-source code used to simulate coupled gas and N-body dynamics in
astrophysical systems, primarily for simulating star cluster formation. Torch combines
different codes via the [AMUSE][1] framework to model gas dynamics, formation of single
and binary stars, stellar evolution and feedback, and collisional N-body dynamics.

See our [website][2] for more information.

The current stable version is [Torch 2.2.3][3].  
The stable versions of Torch live on the [main][4] branch, and ongoing development
happens on the [develop][5] branch.

Quick-start
-----------

The [quickstart][6] is available on our website.

Developing and working with Torch
---------------------------------

Torch is an open-source code with a release brach [main][4] maintained by 
[@torch-sf/teams/developers][17].

Issue reports, enhancements, etc are very welcome.  We will do our best to address
issue reports, pull requests, and inquiries, but due to limited resources we cannot
promise prompt replies.

If you would like to help build Torch please reach out. We have weekly meetings to
discuss progress, future plans on on-going research.

Version philosophy
------------------------
Torch versions are tracked using tags that link to specific commits in the git history.
Releases of Torch follow the version convention torch-vX.Y.Z, where the version X, Y, 
and Z represent major, minor, and hotfix releases, respectively. Feature branches 
should be tagged with appropriate names for use in publications.

AI policy
---------
The Torch developers’ team embraces the AI policies adopted by major peer-reviewed
journals and pre-print services; an example of such policies can be found [here][7] for
the AAS journals. When applicable, individual members of the developers’ team further 
abide by the AI use guidelines of their own institutions. All code in the main and
develop branches is reviewed, tested, and merged directly by members of the
[@torch-sf/teams/developers][17] team. The [@torch-sf/teams/developers][17] team
consists of human researchers with experience running and developing Torch code.

The statement above reflects the workflow of the [@torch-sf/teams/developers][17] and
our policy regarding code in the [main][4] and [develop][5] branches. Although the Torch
code is public and the Torch users’ community is not a formal collaboration, we
recommend that individual users adhere to those principles in their own work.

Acknowledging or citing Torch
-----------------------------
Torch is an open source software and we strongly encourage new users to use and play
with our code. If using Torch or any of our provided resources in work/research
presented in publication, we ask that you cite the relevant papers.

The orginal Torch papers are [Wall et al., 2019][8] and [Wall et al., 2020][9].
Torch also include major upgrades so we ask that you look at and consider whether it is
relevant to cite [Cournoyer-Cloutier et al., 2021][10]; [Polak et al., 2024][11];
[Cournoyer-Cloutier et al, 2025][12]; [Lewis et al., 2025][13]; 
[Appel et al., 2026][14]. In addition, any paper that uses Torch should acknowledge the
FLASH code (as described here: [flash.rochester.edu/site/flashcode.html][15]) and AMUSE
(as described here: [www.amusecode.org/copyright][16]). We also ask that you include
the following statement in the acknowledgements section: 

```
This work made use of the star cluster formation framework Torch \footnote{github.com/torch-sf/torch}.
```

We understand that Torch builds on many different bits and pieces and that it can be
hard to understand exactly what parts you used. Therefore, we encourage new users to
reach out to any of our active developers, attend one of our meetings, and/or join our
slack. We are a friendly group of people who get very happy when people outside the
collaboration want to use our code!

Credits
=======

The Torch code includes contributions by:

* Eric Andersson
* Sabrina Appel
* Claude Cournoyer-Cloutier
* Will Farner
* Joseph Glaser
* Ralf Klessen
* Sean Lewis
* Steve McMillan
* Mordecai-Mark Mac Low
* Andrew Pellegrino
* Brooke Polak
* Simon Portegies Zwart
* Steven Rieder
* Aaron Tran
* Lourens Veen
* Joshua Wall
* Maite Wilhelm

In addition to the FLASH and AMUSE codes, Torch also builds upon software by:

* Christian Bacyznski
* Robi Banerjee
* Moo Kwang Ryan Joung
* Juan Camilo Ibáñez-Mejía
* Daniel Seifried
* Long Wang
* Richard Wünsch

Torch also acknowledges the contributions of everyone who helped build FLASH and AMUSE.

[1]: https://www.amusecode.org/
[2]: https://torch-sf.github.io/
[3]: https://github.com/torch-sf/torch/releases/tag/torch-v2.2.3
[4]: https://github.com/torch-sf/torch/tree/main
[5]: https://github.com/torch-sf/torch/tree/develop
[6]: https://torch-sf.github.io/docs/manual/Quickstart.html
[7]: https://journals.aas.org/news/new-guidelines-on-ai-use-in-aas-journals/
[8]: https://ui.adsabs.harvard.edu/abs/2019ApJ...887...62W/abstract
[9]: https://ui.adsabs.harvard.edu/abs/2020ApJ...904..192W/abstract
[10]: https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.4464C/abstract
[11]: https://ui.adsabs.harvard.edu/abs/2024A%26A...690A..94P/abstract
[12]: https://ui.adsabs.harvard.edu/abs/2025ApJ...990..112C/abstract
[13]: https://ui.adsabs.harvard.edu/abs/2025ApJ...994...69L/abstract
[14]: https://ui.adsabs.harvard.edu/abs/2026AJ....172...12A/abstract
[15]: flash.rochester.edu/site/flashcode.html
[16]: www.amusecode.org/copyright
[17]: https://github.com/orgs/torch-sf/teams/developers
