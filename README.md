# LENT


Description of each folder:

* Figures contains figures organized chronologically by sub-folders "Month_Year" and a Python script that has useful plotting oriented functions "figures_functions.py".

* Magnus_codes contains codes for Magnus evolving a Hamiltonian via the Wegner or relative kinetic energy generator.

* Notes contains pdf notes relevant to some of the tools used in these codes. Add anything that does not qualify as a paper or presentation to this folder. Also contains sub-folders of LaTeX notes (e.g. the draft of our Magnus paper).

* Papers contains all papers relevant to research (i.e. SRG papers, chiral EFT papers, etc.)

* Potentials/vsrg_macos contains all Fortran and Perl codes used to create chiral NN potentials. See the README in this sub-directory for more information on how to generate potentials. Also stores the initial and evolved potentials and momentum arrays and a Python script for loading and saving called "load_save_potentials.py". 

* Presentations contains all the files used for presentations like the annual DNP or NUCLEI meetings.

* SRG_codes contains codes for SRG evolving a Hamiltonian via the Wegner, relative kinetic energy, or block-diagonal generator. Also contains a unitary transformation function within "srg_unitary_transformation.py".


Description of major codes:

* Jupyter notebooks titled with an extension "_figures_v#.ipynb" are used to generate figures which go to the Figures folder.

* "deuteron.py" - Calculates deuteron observables. Work in progress.

* "evolve_hamiltonian.py" - This code SRG or Magnus evolves a specified Hamiltonian using scripts from the SRG_codes or Magnus_codes folder.

* "scratchwork.py" - Script used for testing purposes.
