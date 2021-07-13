 MODULE formfactor

	USE input
	USE math

	IMPLICIT NONE

 CONTAINS


	DOUBLE PRECISION FUNCTION GetFormFactor_1(qMomentum, Isospin)
		INTEGER, INTENT(IN) :: Isospin
		DOUBLE PRECISION, INTENT(IN) :: qMomentum

		INTEGER :: i
		DOUBLE PRECISION :: Argument, BesselFunc, BesselInteger
		DOUBLE PRECISION :: h, res, pi

		DOUBLE PRECISION, ALLOCATABLE :: Fonction(:)

		h = Radius(Npoint) - Radius(Npoint-1)

		pi = 4.0D0*ATAN(1.0D0)

		IF (Isospin == 0) THEN

			ALLOCATE(Fonction(1:Npoint))

			DO i = 1, Npoint

				Argument = qMomentum*Radius(i)
				BesselInteger = j0Spherical(Argument)

				Fonction(i) = BesselInteger * WaveFunVNeut(i) * Radius(i)**2

			END DO
			Fonction(1)=0.0D0

			CALL simps(Fonction, Npoint, h, res)

			DEALLOCATE(Fonction)

			GetFormFactor_1 = 4.0D0*pi*res

		END IF

		IF (Isospin == 1) THEN

			ALLOCATE(Fonction(1:Npoint))

			DO i = 1, Npoint

				Argument = qMomentum*Radius(i)
				BesselInteger = j0Spherical(Argument)

				Fonction(i) = BesselInteger * WaveFunVProt(i) * Radius(i)**2

			END DO
			Fonction(1)=0.0D0

			CALL simps(Fonction, Npoint, h, res)

			DEALLOCATE(Fonction)

			GetFormFactor_1 = 4.0*pi*res

		END IF

		RETURN

	END FUNCTION GetFormFactor_1

	! Subroutine that computes the form factor by taking the
	! Bessel tranform of the nuclear density.
	!
	! Inputs are: - qMomentum
	!             - Isospin
	!

	DOUBLE PRECISION FUNCTION GetFormFactor(qMomentum, Isospin)
		INTEGER, INTENT(IN) :: Isospin
		DOUBLE PRECISION, INTENT(IN) :: qMomentum

		INTEGER :: i
		DOUBLE PRECISION :: Argument, BesselFunc, BesselInteger
		DOUBLE PRECISION :: rk, rip, rkp, h, res, pi

		DOUBLE PRECISION, ALLOCATABLE :: Fonction(:)

		h = Radius(Npoint) - Radius(Npoint-1)

		pi = 4.0D0*ATAN(1.0D0)

		IF (Isospin == 0) THEN

			ALLOCATE(Fonction(1:Npoint))

			DO i = 1, Npoint

				Argument = qMomentum*Radius(i)
				BesselInteger = j0Spherical(Argument)

				Fonction(i) = BesselInteger * DensityNeutron(i) * Radius(i)**2

			END DO

			CALL simps(Fonction, Npoint, h, res)

			DEALLOCATE(Fonction)

			GetFormFactor = 4.0D0*pi*res

		END IF

		IF (Isospin == 1) THEN

			ALLOCATE(Fonction(1:Npoint))

			DO i = 1, Npoint

				Argument = qMomentum*Radius(i)
				BesselInteger = j0Spherical(Argument)

				Fonction(i) = BesselInteger * DensityProton(i) * Radius(i)**2

			END DO

			CALL simps(Fonction, Npoint, h, res)

			DEALLOCATE(Fonction)

			GetFormFactor = 4.0D0*pi*res

		END IF

		IF (Isospin == 2) THEN

			ALLOCATE(Fonction(1:Npoint))

			DO i = 1, Npoint

				Argument = qMomentum*Radius(i)
				BesselInteger = j0Spherical(Argument)

				Fonction(i) = BesselInteger * (DensityNeutron(i) + DensityProton(i)) * Radius(i)**2

			END DO

			CALL simps(Fonction, Npoint, h, res)

			DEALLOCATE(Fonction)

			GetFormFactor = 4.0D0*pi*res

		END IF

		RETURN

	END FUNCTION GetFormFactor

	! Subroutine that computes the form factor by taking the
	! Bessel tranform of the nuclear density.
	!
	! Inputs are: - qMomentum
	!             - Isospin
	!

	DOUBLE PRECISION FUNCTION GetParticleNumber(Isospin)
		INTEGER, INTENT(IN) :: Isospin

		INTEGER :: i
		DOUBLE PRECISION :: h, res, pi

		DOUBLE PRECISION, ALLOCATABLE :: Fonction(:)

		h = Radius(Npoint) - Radius(Npoint-1)

		pi = 4.0D0*ATAN(1.0D0)

		IF (Isospin == 0) THEN

			ALLOCATE(Fonction(1:Npoint))

			DO i = 1, Npoint
				Fonction(i) = DensityNeutron(i) * Radius(i)**2
			END DO

			CALL simps(Fonction, Npoint, h, res)

			DEALLOCATE(Fonction)

			GetParticleNumber = 4.0D0*pi*res

		END IF

		IF (Isospin == 1) THEN

			ALLOCATE(Fonction(1:Npoint))

			DO i = 1, Npoint
				Fonction(i) = DensityProton(i) * Radius(i)**2
			END DO

			CALL simps(Fonction, Npoint, h, res)

			DEALLOCATE(Fonction)

			GetParticleNumber = 4.0D0*pi*res

		END IF

		RETURN

	END FUNCTION GetParticleNumber

	DOUBLE PRECISION FUNCTION GetNorme(Isospin)
		INTEGER, INTENT(IN) :: Isospin

		INTEGER :: i
		DOUBLE PRECISION :: h, res, pi

		DOUBLE PRECISION, ALLOCATABLE :: Fonction(:)

		h = Radius(Npoint) - Radius(Npoint-1)

		pi = 4.0D0*ATAN(1.0D0)

		IF (Isospin == 0) THEN

			ALLOCATE(Fonction(1:Npoint))

			DO i = 1, Npoint
				Fonction(i) = WaveFunVNeut(i)**2 * Radius(i)**2
			END DO

			CALL simps(Fonction, Npoint, h, res)

			DEALLOCATE(Fonction)

			GetNorme = 4.0D0*pi*res

		END IF

		IF (Isospin == 1) THEN

			ALLOCATE(Fonction(1:Npoint))

			DO i = 1, Npoint
				Fonction(i) = WaveFunVProt(i)**2 * Radius(i)**2
			END DO

			CALL simps(Fonction, Npoint, h, res)

			DEALLOCATE(Fonction)

			GetNorme = 4.0D0*pi*res

		END IF

		RETURN

	END FUNCTION GetNorme
	! Subroutine that computes the derivative of the form factor

	DOUBLE PRECISION FUNCTION DerivFormFactor(x, Isospin)
		INTEGER, INTENT(IN) :: Isospin
		DOUBLE PRECISION, INTENT(IN) :: x

		DOUBLE PRECISION :: h = 1.E-3, q1, q2, q3, q4, F1, F2, F3, F4

		q1 = x - 2.0D0*h
		F1 = GetFormFactor(q1, Isospin)

		q2 = x - h
		F2 = GetFormFactor(q2, Isospin)

                q3 = x + h
                F3 = GetFormFactor(q3, Isospin)

                q4 = x + 2.0D0*h
                F4 = GetFormFactor(q4, Isospin)

                DerivFormFactor = (F1 - 8.0D0*F2 + 8.0D0*F3 - F4)/(12.0D0*h)

		RETURN

	END FUNCTION DerivFormFactor

	! Subroutine that computes the actual Helm form factor and returns the
	! surface thickness

	SUBROUTINE FormFactorHelm(qMomentum, qFirstMax, R0, FormHelm, Sigma_2, Isospin)
		INTEGER, INTENT(IN) :: Isospin
		DOUBLE PRECISION, INTENT(IN) :: qMomentum, qFirstMax, R0
		DOUBLE PRECISION, INTENT(OUT) :: FormHelm, Sigma_2

		DOUBLE PRECISION :: Argument, FormFirstMax, j1R0, j1qm, N

                 N = GetParticleNumber(Isospin)

		                             Argument = qMomentum*R0
		IF (ABS(Argument).LT.1.E-14) Argument = 1.D-14

		j1R0 = j1Spherical(Argument)

		                             Argument = qFirstMax*R0
		IF (ABS(Argument).LT.1.E-14) Argument = 1.D-14

		j1qm = j1Spherical(Argument)

		FormFirstMax = GetFormFactor(qFirstMax, Isospin)

		Sigma_2 = (2.0D0/qFirstMax**2) * LOG(ABS(3.0D0*N*j1qm/(qFirstMax*R0*FormFirstMax)))

		FormHelm = 3.0D0/(qMomentum*R0)*j1R0*EXP(-0.5D0*Sigma_2*qMomentum**2)

		RETURN

	END SUBROUTINE FormFactorHelm

	! The first zero of this function will define the diffraction radius for neutrons

	DOUBLE PRECISION FUNCTION Fonction_RootNeut(xArgument, Fderiv)
		DOUBLE PRECISION, INTENT(IN) :: xArgument
		DOUBLE PRECISION, INTENT(OUT) :: Fderiv

		Fderiv = DerivFonction(GetFormFactor, xArgument, 0)
		Fonction_RootNeut = GetFormFactor(xArgument, 0)

		RETURN
	END FUNCTION Fonction_RootNeut

	! The first zero of this function will define the diffraction radius for protons

	DOUBLE PRECISION FUNCTION Fonction_RootProt(xArgument, Fderiv)
		DOUBLE PRECISION, INTENT(IN) :: xArgument
		DOUBLE PRECISION, INTENT(OUT) :: Fderiv

		Fderiv = DerivFonction(GetFormFactor, xArgument, 1)
		Fonction_RootProt = GetFormFactor(xArgument, 1)

		RETURN
	END FUNCTION Fonction_RootProt

	! The first zero of this function will define the position of the first maximum
	! of the form-factor for neutrons

	DOUBLE PRECISION FUNCTION Fonction_RootDerivNeut(xArgument, Fderiv)
		DOUBLE PRECISION, INTENT(IN) :: xArgument
		DOUBLE PRECISION, INTENT(OUT) :: Fderiv

		Fderiv = DerivFonction(DerivFormFactor, xArgument, 0)
		Fonction_RootDerivNeut = DerivFormFactor(xArgument, 0)

		RETURN
	END FUNCTION Fonction_RootDerivNeut

	! The first zero of this function will define the position of the first maximum
	! of the form-factor for protons

	DOUBLE PRECISION FUNCTION Fonction_RootDerivProt(xArgument, Fderiv)
		DOUBLE PRECISION, INTENT(IN) :: xArgument
		DOUBLE PRECISION, INTENT(OUT) :: Fderiv

		Fderiv = DerivFonction(DerivFormFactor, xArgument, 1)
		Fonction_RootDerivProt = DerivFormFactor(xArgument, 1)

		RETURN
	END FUNCTION Fonction_RootDerivProt

	! Generic subroutine that computes the derivative of a function
	! of 2 arguments, the second one being the isospin

	DOUBLE PRECISION FUNCTION DerivFonction(fonctionName, x, Isospin)
		DOUBLE PRECISION, EXTERNAL :: fonctionName
		INTEGER, INTENT(IN) :: Isospin
		DOUBLE PRECISION, INTENT(IN) :: x
		DOUBLE PRECISION :: h = 1.E-3, q1, q2, q3, q4, F1, F2, F3, F4

		q1 = x - 2.0D0*h
		F1 = fonctionName(q1, Isospin)

		q2 = x - h
		F2 = fonctionName(q2, Isospin)

                q3 = x + h
                F3 = fonctionName(q3, Isospin)

                q4 = x + 2.0D0*h
                F4 = fonctionName(q4, Isospin)

		DerivFonction = (F1 - 8.0D0*F2 + 8.0D0*F3 - F4)/(12.0D0*h)

	RETURN
	END FUNCTION DerivFonction

	! Spherical Bessel function at n=0
	DOUBLE PRECISION FUNCTION j0Spherical(x)
		DOUBLE PRECISION, INTENT(IN) :: x

		j0Spherical = SIN(x)/x

		RETURN
	END FUNCTION j0Spherical

	! Spherical Bessel function at n=1
	DOUBLE PRECISION FUNCTION j1Spherical(x)
		DOUBLE PRECISION, INTENT(IN) :: x

		j1Spherical = SIN(x)/x**2 - COS(x)/x

		RETURN
	END FUNCTION j1Spherical


END MODULE formfactor
