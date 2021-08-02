Here we describe modifications and tests to have only one pion exchange.

If we use Evgeny's potentials, then we can use LO and turn off the S-wave
constants.
  * in ichiral_ope.f90, we modified LO to put 0 in front of the constants.
  * we also introduced ostat=-1 to allow for unregulated and nonlocal regulator.

  kvnn = 35 --> unregulated OPE
  kvnn = 36 non-local Fermi-Dirac regulator with \Lambda = 500 MeV and \eps_FD = 0.5 fm^{-1}
  kvnn = 37,38,39 non-local regulator with e^{-((p^2+p'^2)/\Lambda^2)^nexp}
          with nexp = 1,2,3


The diagonal matrix elements for 3D2, 3F3, and 3G4 are consistent with plots
 in a talk by Evgeny.

For plot comparisons:
  * Look up the Elab to p conversion for Evgeny.

  * We use h^2/M = 41.47105 MeV-fm^2

  * To convert potential from fm to GeV^{-2}, use
        41.47105 * 1/197.33 * 1/(.19733)^2 = 5.39716
   * To convert k_rel to Elab (nonrel for now):  Elab = k_rel^2 * 82.9421
