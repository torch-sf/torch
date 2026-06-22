.. Torch-SF documentation master file, created by
   sphinx-quickstart on Thu Oct 23 16:58:49 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Torch-SF
========

This website is still under construction. The quickstart notes are up to date.

Torch is an open-source code used to simulate coupled gas and N-body dynamics 
in astrophysical systems, primarily for simulating star cluster formation. 
Torch combines different codes via the `AMUSE`_ framework to model gas 
dynamics, formation of single and binary stars, stellar evolution and feedback, 
and collisional N-body dynamics.

Stable versions of torch live on the `main`_ branch. Ongoing development 
happens on the `develop`_ branch.

.. _AMUSE: https://www.amusecode.org/
.. _main: https://github.com/torch-sf/torch
.. _develop: https://github.com/torch-sf/torch/tree/develop

.. toctree::
   :maxdepth: 2
   :caption: User guide
   
   manual/Quickstart.md
   manual/Installation.md
   manual/FAQ.md
.. manual/Torch_Parameters.md
   manual/Flash_Parameters.md
   manual/Known_Systems.md
   manual/Known_Issues.md

.. toctree::
   :maxdepth: 2
   :caption: Wiki
   
   wiki/Torch.md
.. wiki/Magnetohydrodynamics.md
   wiki/Cooling_Heating.md
   wiki/Star_Formation.md
   wiki/Stellar_Feedback.md
   wiki/Radiation.md
   wiki/Nbody.md
   wiki/Initial_Conditions.md

.. toctree::
   :maxdepth: 2
   :caption: Developer documentation

   developer/Developer_Guide.md
.. developer/Amuse_Guide.md
   developer/Flash_Guide.md
   developer/In_Progress.md
   developer/modules

