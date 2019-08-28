# LENT

This repository contains codes and files relevant to SRG work done by A.T.


Description of each folder:

* __Figures__ contains research/results figures and a Python script figures_functions.py that has useful plotting oriented functions.

* __Magnus_codes__ contains codes for Magnus-evolving a Hamiltonian via the Wegner or relative kinetic energy generator. Also contains a script for evolving via the "split thing" approach. See magnus_split.py for details.

* __Notes__ contains notes relevant to some of the tools used in these codes. Add anything that does not qualify as a paper or presentation to this folder. Also contains sub-folders of LaTeX notes (e.g. the draft of the Magnus paper, candidacy paper, etc.)

* __Old_codes__ stores all out-dated codes. See the README stored in this folder for documentation.

* __Papers__ contains all relevant papers.

* __Potentials/vsrg_macos__ contains all Fortran and Perl codes used to create chiral NN potentials. See the README in this sub-directory for more information on how to generate potentials. Also stores the initial and evolved potentials, momentum and weights arrays, and a Python script for loading and saving potentials called load_save_potentials.py. 

* __Presentations__ contains all the files used for presentations for conferences or meetings.

* __SRG_codes__ contains codes for SRG evolving a Hamiltonian via the Wegner, relative kinetic energy, or block-diagonal generator. Also contains a unitary transformation function within srg_unitary_transformation.py.


Description of major codes:

* Jupyter notebooks titled with an extension fig.ipynb are used to generate figures which go to the __Figures__ folder.

* phase_shifts.py - Calculates phase shifts given an NN potential. Work in progress.

* observables.py - Calculates NN observables. Work in progress.

* operators.py - Creates momentum-space operators for NN observables. Work in progress.

* run_magnus.py - Magnus-evolves a specified Hamiltonian using scripts from the __Magnus_codes__ folder.

* run_srg.py - SRG-evolves a specified Hamiltonian using scripts from the __SRG_codes__ folder.

* test.py - Script used for testing purposes.
