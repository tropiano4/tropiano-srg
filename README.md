# LENT

This repository contains codes and files relevant to SRG work done by A.T.


Description of each folder:

* __Densities__ contains codes for computing the nucleonic densities of different nuclei.

* __Figures__ contains research/results figures and a Python script figures_functions.py that has useful plotting oriented functions.

* __Magnus__ contains codes for SRG-evolving a Hamiltonian via the Wegner or relative kinetic energy generator and using the Magnus expansion.

* __Notes__ contains notes relevant to some of the tools used in these codes. Add anything that does not qualify as a paper or presentation to this folder. Also contains sub-folders of LaTeX notes (e.g. the arXiv version of the operator evolution paper, candidacy paper, etc.) and notes taken on iPad.

* __Old_codes__ stores all out-dated codes.

* __Papers__ contains all relevant papers. These papers are organized by the BibTeX article name (e.g., Tropiano:2018quk would be named Tropiano_2018quk.pdf). 

* __Potentials__ contains all Fortran and Perl codes used to create chiral NN potentials. See the README for more information on how to generate potentials. Also stores the initial and evolved potentials, momenta and weights, and a Python script for loading and saving potentials called vnn.py. 

* __Presentations__ contains all the files used for presentations for conferences or meetings.

* __SRG__ contains codes for SRG evolving a Hamiltonian via the Wegner, relative kinetic energy, or block-diagonal generator. Also contains a U(k,k') transformation function within srg_unitary_transformation.py.


Description of major codes:

* Jupyter notebooks titled with an extension fig.ipynb are used to generate figures which go to the __Figures__ folder.

* densities.py - Loads nucleonic densities from the __Densities__ folder. These densities are used in conjunction with snmd.py, pmd.py, and dmd.py.

* observables.py - Calculates NN observables and wave functions.

* operators.py - Creates momentum-space operators.

* run_magnus.py - Magnus-evolves a specified Hamiltonian using scripts from the __Magnus__ folder.

* run_srg.py - SRG-evolves a specified Hamiltonian using scripts from the __SRG__ folder.

* snmd.py, pmd.py, dmd.py - Calculates single-nucleon, pair, and deuteron momentum distributions in LDA.

* test_script.py - Script used for testing purposes.
