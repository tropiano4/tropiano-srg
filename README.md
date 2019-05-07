# LENT


Description of each folder:

* Figures contains figures organized chronologically by sub-folders "Month_Year". This folder still needs work for better organization. Can add README to each sub-folder to explain each of the figures if the title of the figure is not sufficient.

* Fortran_test_codes contains Fortran codes used to test SRG or Magnus algorithms.

* Notes contains pdf notes relevant to some of the tools used in these codes. Add anything that does not qualify as a paper or presentation to this folder. Also contains sub-folders of LaTeX notes (e.g. the draft of our Magnus paper).

* Papers contains all papers relevant to research (i.e. SRG papers, chiral EFT papers, etc.)

* Potentials/vsrg_macos contains all Fortran and Perl codes used to create chiral NN potentials. See the README in this sub-directory for more information on how to generate potentials. Also stores the initial and evolved potentials and momentum arrays and a Python script for loading and saving called "load_save_potentials.py". 

* Presentations contains all the files used for presentations like the annual DNP or NUCLEI meetings.

Description of major codes:

* Jupyter notebooks titled with an extension "_figures_v#.ipynb" are used to generate figures which go to the Figures folder.

* evolve_hamiltonian.py - This code SRG or Magnus evolves a specified Hamiltonian using scripts from the SRG_codes or Magnus_codes folder.

* scratchwork.py - Script used for testing purposes.
