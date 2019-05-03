# LENT

This repository is for research done during my (Anthony Tropiano) PhD at OSU.

Description of each folder:

* Figures contains figures organized chronologically by sub-folders "Month_Year". This folder still needs work for better organization. Can add README to each sub-folder to explain each of the figures if the title of the figure is not sufficient.

* Notes contains pdf notes relevant to some of the tools used in these codes. Add anything that does not qualify as a paper or presentation to this folder. Also contains sub-folders of LaTeX notes (e.g. the draft of our Magnus paper).

* Papers contains all papers relevant to research (i.e. SRG papers, chiral EFT papers, etc.)

* Potentials/vsrg_macos contains all Fortran and Perl codes used to create chiral NN potentials. Also stores the initial and evolved potentials and momentum arrays. See the README in this sub-directory for more information on how to generate potentials.

* Presentations contains all the files used for presentations like the annual DNP or NUCLEI meetings.

Description of major codes:

* Jupyter notebooks titled "NAME_figures_v#.ipynb" are used to generate figures which go to the Figures folder.

* evolve_hamiltonian.py - This code SRG or Magnus evolves a specified Hamiltonian using scripts from the SRG_codes or Magnus_codes folder.

* scratchwork.py - Script used for testing purposes.
