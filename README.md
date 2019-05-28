# LENT

This repository contains codes and files relevant to SRG work done by A.T.


Description of each folder:

* Figures contains research/results figures organized chronologically by sub-folders "Month_Year" and a Python script "figures_functions.py" that has useful plotting oriented functions.

* Magnus_codes contains codes for Magnus evolving a Hamiltonian via the Wegner or relative kinetic energy generator.

* Notes contains pdf notes relevant to some of the tools used in these codes. Add anything that does not qualify as a paper or presentation to this folder. Also contains sub-folders of LaTeX notes (e.g. the draft of our Magnus paper).

* Papers contains all relevant papers and a .bib file for LaTeX citing purposes.

* Potentials/vsrg_macos contains all Fortran and Perl codes used to create chiral NN potentials. See the README in this sub-directory for more information on how to generate potentials. Also stores the initial and evolved potentials, momentum and weights arrays, and a Python script for loading and saving called "load_save_potentials.py". 

* Presentations contains all the files used for presentations for conferences or meetings.

* SRG_codes contains codes for SRG evolving a Hamiltonian via the Wegner, relative kinetic energy, or block-diagonal generator. Also contains a unitary transformation function within "srg_unitary_transformation.py".


Description of major codes:

* Jupyter notebooks titled with an extension "_figures.ipynb" are used to generate figures which go to the Figures folder.

* "deuteron.py" - Calculates deuteron observables. Work in progress.

* "run_magnus.py" - This code Magnus evolves a specified Hamiltonian using scripts from the Magnus_codes folder. Work in progress.

* "run_srg.py" - This code SRG evolves a specified Hamiltonian using scripts from the SRG_codes folder.

* "phase_shifts.py" - Calculates phase shifts given an NN potential. Work in progress.

* "test.py" - Script used for testing purposes.
