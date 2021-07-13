PROGRAM main

	USE input
	USE formfactor
	USE math

	IMPLICIT NONE

	LOGICAL :: do_helm = .False.

	INTEGER :: cycles_in, cycles_out, cycles_rate, i
	INTEGER :: Nmax = 200, Nstep = 0

	DOUBLE PRECISION :: fatol, x0, f0, norme_n, norme_p, q_max
	DOUBLE PRECISION :: qMomentum, qFirstMax , FormHelm, R0Diffraction, Sigma_2
	DOUBLE PRECISION :: R0Diffraction_n, Sigma_2_n, R0Diffraction_p, Sigma_2_p
	DOUBLE PRECISION :: Rhelm_n, Rhelm_p

	DOUBLE PRECISION, ALLOCATABLE :: FormFac(:), Qmoment(:)

	CHARACTER ( LEN = 10 ) namvar

        !---------------------------------------------------------------------!
        !    								      !
        !                 BEGINNING OF THE PROGRAM                            !
        !    								      !
        !---------------------------------------------------------------------!

	CALL SYSTEM_CLOCK(cycles_in)

	CALL ReadDensity()
	!CALL ReadWaveFunction()
	!CALL PrintOut(Radius, DensityNeutron)

	ALLOCATE(FormFac(1:Npoint))
	ALLOCATE(Qmoment(1:Npoint))

	q_max = 5.0D0
	DO i = 1, Npoint

		qMomentum = REAL(i,KIND=8)/REAL(Npoint,KIND=8)*q_max

		Qmoment(i) = qMomentum
		FormFac(i) = GetFormFactor(qMomentum, 2)

	END DO

	CALL PrintOut(Qmoment, FormFac)

	norme_n = GetNorme(0)
	norme_p = GetNorme(1)
	write(6,'("Norme neutrons = ",f20.14)') norme_n
	write(6,'("Norme protons = ",f20.14)') norme_p

	IF (do_helm) THEN

		fatol = 1.E-6

		! Finds the first zero of the Helm form factor function (NEUTRONS)

		x0 = 0.5

		CALL Newton(Fonction_RootNeut, Nmax, Nstep, fatol, x0, f0, 0)

		write(*,'("First zero: ",f12.8)') x0

		R0Diffraction = 4.49341/x0

		! Find the first maximum of the Helm form factor by looking at the first
		! zero of its derivative (NEUTRONS)

		x0 = x0 + 0.02

		CALL Newton(Fonction_RootDerivNeut, Nmax, Nstep, fatol, x0, f0, 1)

		write(*,'("First Maximum: ",f12.8)') x0

		! Computes the diffraction radius R0Diffraction and surface thickness
		! Sigma_2 (NEUTRONS)

		qFirstMax = x0
		CALL FormFactorHelm(qMomentum, qFirstMax, R0Diffraction, FormHelm, Sigma_2, 0)

		R0Diffraction_n = R0Diffraction
		Sigma_2_n = Sigma_2
        	Rhelm_n = SQRT(R0Diffraction**2+5.0*Sigma_2)

		! Finds the first zero of the Helm form factor function (PROTONS)

		x0 = 0.5

		CALL Newton(Fonction_RootProt, Nmax, Nstep, fatol, x0, f0, 0)

        	write(*,'("First zero: ",f12.8)') x0

		R0Diffraction = 4.49341/x0

		! Find the first maximum of the Helm form factor by looking at the first
		! zero of its derivative (PROTONS)

		x0 = x0 + 0.1

		CALL Newton(Fonction_RootDerivProt, Nmax, Nstep, fatol, x0, f0, 1)

		write(*,'("First maximum: ",f12.8)') x0

		! Computes the diffraction radius R0Diffraction and surface thickness
		! Sigma_2 (PROTONS)

		qFirstMax = x0
		CALL FormFactorHelm(qMomentum, qFirstMax, R0Diffraction, FormHelm, Sigma_2, 1)

		R0Diffraction_p = R0Diffraction
		Sigma_2_p = Sigma_2
		Rhelm_p = SQRT(R0Diffraction**2+5.0*Sigma_2)

		WRITE(*,'("Neutrons:",3E24.16)') Rhelm_n, R0Diffraction_n, SQRT(Sigma_2_n)
		WRITE(*,'("Protons: ",3E24.16)') Rhelm_p, R0Diffraction_p, SQRT(Sigma_2_p)

	END IF

	DEALLOCATE(FormFac)
	DEALLOCATE(Qmoment)

	CALL SYSTEM_CLOCK(cycles_out, cycles_rate)

	STOP

END PROGRAM main
