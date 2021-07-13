# Generating densities with Gogny code from Nicolas Schunck

OK, so a few words of explanations:

* There is only one branch;
* You should look at what is under Bessel/, this is the ‘vanilla’ spherical Gogny code;
* The Makefile should be rather straightforward, it’s preset for gfortran and should work out of the box since there is no external library;
* After compilation, just edit the file data/input.txt. Apart from the obvious like Z, N, b0 (oscillator length), N0 (number of shells), you may only want to change the last 4 lines that activate/deactivate specific terms of the potential in the calculation of the mean field or pairing field, and of the expectation value. Value 0 means the corresponding term is not included:
  * Coulomb: 1 = direct term only, 2 = direct+exchange, 3=direct+exchange+pairing
  * CenterofMass: this is the 2-body center of mass: 1 = mean field contribution, 2 = mean-field + pairing
  * SpinOrbit: 1 = mean-field, 2 = mean-field + pairing
* Run the code simply by typing: main > main.out
* Some additional output is written in data/, the densities are in DensityQP.dat (r, rho_p, rho_n). These are the usual radial densities.
* It is safer to remove the matrix elements of the potentials, i.e., all the files *HO.txt before running a new calculation.

Crank up the number of oscillator shells. N0 = 16 seems good.

Use b = np.sqrt(0.90*A**(1/3) + 0.7) fm for the oscillator parameter.
  (This comes from some old paper.)
 

