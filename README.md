# LENT

This repository contains codes and files relevant to SRG work done by A.T.


Description of each folder:

* __Figures__ contains research/results figures and a Python script _figures_functions.py_ that has useful plotting oriented functions.

* __Magnus_codes__ contains codes for Magnus-evolving a Hamiltonian via the Wegner or relative kinetic energy generator. Also contains a script for evolving via the "split thing" approach. See _magnus_split.py_ for details.

* __Notes__ contains notes relevant to some of the tools used in these codes. Add anything that does not qualify as a paper or presentation to this folder. Also contains sub-folders of LaTeX notes (e.g. the draft of the Magnus paper, candidacy paper, etc.)

* __Old_codes__ stores all out-dated codes. See the README stored in this folder for documentation.

* __Papers__ contains all relevant papers.

* __Potentials/vsrg_macos__ contains all Fortran and Perl codes used to create chiral NN potentials. See the README in this sub-directory for more information on how to generate potentials. Also stores the initial and evolved potentials, momentum and weights arrays, and a Python script for loading and saving potentials called _load_save_potentials.py_. 

* __Presentations__ contains all the files used for presentations for conferences or meetings.

* __SRG_codes__ contains codes for SRG evolving a Hamiltonian via the Wegner, relative kinetic energy, or block-diagonal generator. Also contains a unitary transformation function within _srg_unitary_transformation.py_.


Description of major codes:

* Jupyter notebooks titled with an extension _figures.ipynb_ are used to generate figures which go to the __Figures__ folder.

* _phase_shifts.py_ - Calculates phase shifts given an NN potential. Work in progress.

* _observables.py_ - Calculates NN observables. Work in progress.

* _operators.py_ - Creates momentum-space operators for NN observables. Work in progress.

* _run_magnus.py_ - Magnus-evolves a specified Hamiltonian using scripts from the __Magnus_codes__ folder.

* _run_srg.py_ - SRG-evolves a specified Hamiltonian using scripts from the __SRG_codes__ folder.

* _test.py_ - Script used for testing purposes.
