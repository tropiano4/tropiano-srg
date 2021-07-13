 !----------------------------------------------------------------------------------------------!
 !    This module defines an object called "selfconsistencymethodproj" and contains a set of	!
 !    functions that calculate various energies (brink-boeker, spin-orbit, coulomb, etc.).	!
 !    It also contains a few utilities to display the status of a fiven hfb iteration on screen.!
 !----------------------------------------------------------------------------------------------!

 MODULE selfc_proj

	USE input
	USE nucleus
	USE symden_proj
	USE symke2b
	USE symvbb
	USE symvc
	USE symvls
	USE symgdd_proj

	IMPLICIT NONE

	! This type is used to contain at once all matrix elements as well as the density

	TYPE SelfConsistencyMethodProj
		TYPE (SymKineticEnergy2Body) vEkCMph
		TYPE (SymEk2pp) vEkCMpp
		TYPE (SymVBBph) vBBph
		TYPE (SymVBBpp) vBBpp
		TYPE (SymVCph) vCph
		TYPE (SymVCpp) vCpp
		TYPE (SymVLSph) vLSph
		TYPE (SymVLSpp) vLSpp
		TYPE (SymGDDph) gDDph
		TYPE (SymDensityProj), POINTER :: density
	END TYPE

 CONTAINS

	SUBROUTINE SelfConsistencyMethodProj_new(consistency, density)
		TYPE (SelfConsistencyMethodProj), INTENT(INOUT) :: consistency
		TYPE (SymDensityProj), TARGET, INTENT(IN) :: density

		DOUBLE PRECISION b

		b = Nucleus_get_b(density%nucleus)

		CALL SymKineticEnergy2Body_new(consistency%vEkCMph, consistency%vEkCMpp)
		CALL SymVBBph_new(consistency%vBBph, b)
		CALL SymVBBpp_new(consistency%vBBpp, b)
		CALL SymVCph_new(consistency%vCph)
		CALL SymVCpp_new(consistency%vCpp)
		CALL SymVLSph_new(consistency%vLSph)
		CALL SymVLSpp_new(consistency%vLSpp)
		CALL SymGDDph_new(consistency%gDDph)

		consistency%density => density

		RETURN
	END SUBROUTINE SelfConsistencyMethodProj_new

	SUBROUTINE SelfConsistencyMethodProj_store_eHFB(consistency)
		TYPE (SelfConsistencyMethodProj), INTENT(INOUT) :: consistency

		DOUBLE PRECISION eHF
		DOUBLE PRECISION kinetic_energy, kinetic_CM_energy
		DOUBLE PRECISION local_energy_BB, exchange_energy_BB ! Brink-Booker
		DOUBLE PRECISION local_energy_Coulomb, exchange_energy_Coulomb
		DOUBLE PRECISION dd_energy ! Density Dependent
		DOUBLE PRECISION ls_energy

		DOUBLE PRECISION pairing
		DOUBLE PRECISION pairing_Ek
		DOUBLE PRECISION pairing_BB ! BrinkBooker
		DOUBLE PRECISION pairing_Coulomb
		DOUBLE PRECISION pairing_LS !

		kinetic_energy          = SelfConsistencyMethodProj_get_Ek(consistency)
		kinetic_CM_energy       = SelfConsistencyMethodProj_get_EkCM(consistency)
		local_energy_BB         = SelfConsistencyMethodProj_get_LocalBBEnergy(consistency)
 		local_energy_Coulomb    = SelfConsistencyMethodProj_get_LocalCoulombEnergy(consistency)
		exchange_energy_BB      = SelfConsistencyMethodProj_get_ExchangeBBEnergy(consistency)
		exchange_energy_Coulomb = SelfConsistencyMethodProj_get_ExchangeCoulombEnergy(consistency)
		dd_energy               = SelfConsistencyMethodProj_get_DDEnergy(consistency)
		ls_energy               = SelfConsistencyMethodProj_get_LSEnergy(consistency)

		eHF = kinetic_energy + kinetic_CM_energy + local_energy_BB + exchange_energy_BB &
                + local_energy_Coulomb + exchange_energy_Coulomb + dd_energy + ls_energy

		pairing_Ek      = SelfConsistencyMethodProj_get_EkPairing(consistency)
		pairing_BB      = SelfConsistencyMethodProj_get_BBPairing(consistency)
		pairing_Coulomb = SelfConsistencyMethodProj_get_CoulombPairing(consistency)
		pairing_LS      = SelfConsistencyMethodProj_get_LSPairing(consistency)

		pairing = pairing_BB + pairing_Coulomb + pairing_LS + pairing_Ek

		consistency%density%nucleus%eHFB = eHF + pairing
		consistency%density%nucleus%pairing = pairing

		RETURN
	END SUBROUTINE SelfConsistencyMethodProj_store_eHFB

	! One-body Kinetic Energy
	FUNCTION SelfConsistencyMethodProj_get_Ek(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_Ek
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		INTEGER :: A, ta
		DOUBLE PRECISION :: factor, b
		DOUBLE PRECISION, DIMENSION(0:1) :: Ek

		A = Nucleus_get_A(consistency%density%nucleus)
                b = Nucleus_get_b(consistency%density%nucleus)

                factor = (1.0D0 - (1.0D0 / DBLE(A)))

		DO ta = 0, 1
			Ek(ta) = factor * (EkField * consistency%density%field%rho%p(ta))
		END DO

		SelfConsistencyMethodProj_get_Ek = Ek(0) + Ek(1)

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_Ek

	! Two-body Kinetic Energy
	FUNCTION SelfConsistencyMethodProj_get_EkCM(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_EkCM
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		INTEGER :: A, ProjectionOn = 0
		TYPE (SymHartreeFockFieldProj) :: HF_Gamma

		CALL SymHartreeFockFieldProj_new(HF_Gamma)

		! Protons
		CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma, consistency%vEkCMph, consistency%density%field%rho, 0, 0, ProjectionOn)
		! Neutrons
		CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma, consistency%vEkCMph, consistency%density%field%rho, 1, 1, ProjectionOn)

		A = Nucleus_get_A(consistency%density%nucleus)

		IF (A .LE. 1) STOP "Abortado"

		IF (switch_CM .GE. 1) THEN
			SelfConsistencyMethodProj_get_EkCM = (0.5D0 / DBLE(A)) * (HF_Gamma * consistency%density%field%rho)
		ELSE
			SelfConsistencyMethodProj_get_EkCM = 0.0D0
		END IF

		CALL SymHartreeFockFieldProj_del(HF_Gamma)
		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_EkCM

	! Energia local de Brink-Booker
	FUNCTION SelfConsistencyMethodProj_get_LocalBBEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_LocalBBEnergy
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		TYPE (SymHartreeFockFieldProj) :: HF_Gamma
		INTEGER :: ProjectionOn = 0

		CALL SymHartreeFockFieldProj_new(HF_Gamma)

		! Protons
		CALL SymVBBph_get_LocalGamma(HF_Gamma, consistency%vBBph, consistency%density%field%rho, 0, 0, ProjectionOn)
		! Neutrons
		CALL SymVBBph_get_LocalGamma(HF_Gamma, consistency%vBBph, consistency%density%field%rho, 1, 1, ProjectionOn)

		SelfConsistencyMethodProj_get_LocalBBEnergy = 0.5D0 * (HF_Gamma * consistency%density%field%rho)

		CALL SymHartreeFockFieldProj_del(HF_Gamma)

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_LocalBBEnergy

	! Energia de intercambio de Brink-Booker
	FUNCTION SelfConsistencyMethodProj_get_ExchangeBBEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_ExchangeBBEnergy
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		TYPE (SymHartreeFockFieldProj) :: HF_Gamma
		INTEGER :: ProjectionOn = 0

		CALL SymHartreeFockFieldProj_new(HF_Gamma)

		! Protons
		CALL SymVBBph_get_ExchangeGamma(HF_Gamma, consistency%vBBph, consistency%density%field%rho, 0, 0, ProjectionOn)
		! Neutrons
		CALL SymVBBph_get_ExchangeGamma(HF_Gamma, consistency%vBBph, consistency%density%field%rho, 1, 1, ProjectionOn)

		SelfConsistencyMethodProj_get_ExchangeBBEnergy = 0.5D0 * (HF_Gamma * consistency%density%field%rho)

		CALL SymHartreeFockFieldProj_del(HF_Gamma)

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_ExchangeBBEnergy

	! Coulomb local energy
	FUNCTION SelfConsistencyMethodProj_get_LocalCoulombEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_LocalCoulombEnergy
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		TYPE (SymHartreeFockFieldProj) :: HF_Gamma
		INTEGER :: ProjectionOn = 0

		CALL SymHartreeFockFieldProj_new(HF_Gamma)

		! Protons
		CALL SymVCph_get_LocalGamma(HF_Gamma, consistency%vCph, consistency%density%field%rho, 0, 0, ProjectionOn)
		! Neutrons
		CALL SymVCph_get_LocalGamma(HF_Gamma, consistency%vCph, consistency%density%field%rho, 1, 1, ProjectionOn)

		IF (switch_Coulomb .GE. 1) THEN
			SelfConsistencyMethodProj_get_LocalCoulombEnergy = 0.5D0 * (HF_Gamma * consistency%density%field%rho)
		ELSE
			SelfConsistencyMethodProj_get_LocalCoulombEnergy = 0.0D0
		END IF

		CALL SymHartreeFockFieldProj_del(HF_Gamma)

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_LocalCoulombEnergy

	! Exchange Coulomb Energy
	FUNCTION SelfConsistencyMethodProj_get_ExchangeCoulombEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_ExchangeCoulombEnergy
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		TYPE (SymHartreeFockFieldProj) :: HF_Gamma
		INTEGER :: ProjectionOn = 0

		CALL SymHartreeFockFieldProj_new(HF_Gamma)

		! Protons
		CALL SymVCph_get_ExchangeGamma(HF_Gamma, consistency%vCph, consistency%density%field%rho, 0, 0, ProjectionOn)
		! Neutrons
		CALL SymVCph_get_ExchangeGamma(HF_Gamma, consistency%vCph, consistency%density%field%rho, 1, 1, ProjectionOn)

		IF (switch_Coulomb .GE. 2) THEN
			SelfConsistencyMethodProj_get_ExchangeCoulombEnergy = 0.5D0 * (HF_Gamma * consistency%density%field%rho)
		ELSE
			SelfConsistencyMethodProj_get_ExchangeCoulombEnergy = 0.0D0
		END IF

		CALL SymHartreeFockFieldProj_del(HF_Gamma)

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_ExchangeCoulombEnergy

	! Energia dependiente de la densidad
	FUNCTION SelfConsistencyMethodProj_get_DDEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_DDEnergy
		TYPE (SelfConsistencyMethodProj), INTENT(INOUT) :: consistency

		DOUBLE PRECISION :: Gauge

		Gauge = consistency%density%field%rho%GaugeAngle(0)

		IF (Basis .EQ. 1) THEN
			CALL SymGDDph_make_DD(consistency%gDDph, consistency%density%field%rho, Gauge)
		ELSE
			CALL Make_DenGenFun(consistency%gDDph, consistency%density%field%rho, Gauge)
		END IF

		IF (switch_DD .GE. 1) THEN
			SelfConsistencyMethodProj_get_DDEnergy = SymGDDph_get_edd(consistency%gDDph)
		ELSE
			SelfConsistencyMethodProj_get_DDEnergy = 0.0D0
		END IF

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_DDEnergy

	! Spin-Orbit Energy
	FUNCTION SelfConsistencyMethodProj_get_LSEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_LSEnergy
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		TYPE (SymHartreeFockFieldProj) :: HF_Gamma
		INTEGER :: ProjectionOn = 0

		CALL SymHartreeFockFieldProj_new(HF_Gamma)

		! Protons
		CALL SymVLSph_get_Gamma(HF_Gamma, consistency%vLSph, consistency%density%field%rho, 0, 0, ProjectionOn)
		! Neutrons
		CALL SymVLSph_get_Gamma(HF_Gamma, consistency%vLSph, consistency%density%field%rho, 1, 1, ProjectionOn)

		IF (switch_LS .GE. 1) THEN
			SelfConsistencyMethodProj_get_LSEnergy = 0.5D0* (HF_Gamma * consistency%density%field%rho)
		ELSE
			SelfConsistencyMethodProj_get_LSEnergy = 0.0D0
		END IF

		CALL SymHartreeFockFieldProj_del(HF_Gamma)

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_LSEnergy

	! Center of Mass pairing
	FUNCTION SelfConsistencyMethodProj_get_EkPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_EkPairing
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

                INTEGER :: A
		DOUBLE PRECISION :: factor

		TYPE (SymHartreeFockFieldProj) :: delta

		CALL SymHartreeFockFieldProj_new(delta)

                A = Nucleus_get_A(consistency%density%nucleus)

		! Protons
		CALL SymKineticEnergy2Body_get_Delta(delta, consistency%vEkCMpp, consistency%density%field%kap, 0, 0, 0)
		! Neutrons
		CALL SymKineticEnergy2Body_get_Delta(delta, consistency%vEkCMpp, consistency%density%field%kap, 1, 1, 1)

                factor = (1.0D0 / A)

		IF (switch_CM .EQ. 1) THEN
			SelfConsistencyMethodProj_get_EkPairing = 0.5D0 * factor * (consistency%density%field%kap * delta)
		ELSE
			SelfConsistencyMethodProj_get_EkPairing = 0.0D0
		END IF

		CALL SymHartreeFockFieldProj_del(delta)

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_EkPairing

	! Brink-Boeker pairing
	FUNCTION SelfConsistencyMethodProj_get_BBPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_BBPairing
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		TYPE (SymHartreeFockFieldProj) :: delta

		CALL SymHartreeFockFieldProj_new(delta)

		! Protons
		CALL SymVBBpp_get_Delta(delta, consistency%vBBpp, consistency%density%field%kap, 0, 0)
		! Neutrons
		CALL SymVBBpp_get_Delta(delta, consistency%vBBpp, consistency%density%field%kap, 1, 1)

		SelfConsistencyMethodProj_get_BBPairing = 0.5D0 * (consistency%density%field%kap * delta)

		CALL SymHartreeFockFieldProj_del(delta)

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_BBPairing

	! Coulomb pairing
	FUNCTION SelfConsistencyMethodProj_get_CoulombPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_CoulombPairing
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		TYPE (SymHartreeFockFieldProj) :: delta

		CALL SymHartreeFockFieldProj_new(delta)

		! Protons
		CALL SymVCpp_get_Delta(delta, consistency%vCpp, consistency%density%field%kap, 0, 0)
		! Neutrons
		CALL SymVCpp_get_Delta(delta, consistency%vCpp, consistency%density%field%kap, 1, 1)

		IF (switch_LS .GE. 1) THEN
			SelfConsistencyMethodProj_get_CoulombPairing = 0.5D0 * (consistency%density%field%kap * delta)
		ELSE
			SelfConsistencyMethodProj_get_CoulombPairing = 0.0D0
		END IF

		CALL SymHartreeFockFieldProj_del(delta)

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_CoulombPairing

	! Spin-orbit pairing
	FUNCTION SelfConsistencyMethodProj_get_LSPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_LSPairing
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		TYPE (SymHartreeFockFieldProj) :: delta

		CALL SymHartreeFockFieldProj_new(delta)

		! Protons
		CALL SymVLSpp_get_Delta(delta, consistency%vLSpp, consistency%density%field%kap, 0, 0)
		! Neutrons
		CALL SymVLSpp_get_Delta(delta, consistency%vLSpp, consistency%density%field%kap, 1, 1)

		IF (switch_LS .GE. 1) THEN
			SelfConsistencyMethodProj_get_LSPairing = 0.5D0 *(consistency%density%field%kap * delta)
		ELSE
			SelfConsistencyMethodProj_get_LSPairing = 0.0D0
		END IF

		CALL SymHartreeFockFieldProj_del(delta)

		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_LSPairing

	FUNCTION SelfConsistencyMethodProj_accuracy(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_accuracy
		TYPE (SelfConsistencyMethodProj), INTENT(INOUT) :: consistency

		DOUBLE PRECISION old_eHFB

		old_eHFB = consistency%density%nucleus%eHFB

		CALL SelfConsistencyMethodProj_store_eHFB(consistency)

		SelfConsistencyMethodProj_accuracy = ABS(old_eHFB - consistency%density%nucleus%eHFB)
		RETURN
	END FUNCTION SelfConsistencyMethodProj_accuracy

	SUBROUTINE SelfConsistencyMethodProj_show_Status(consistency)
		TYPE (SelfConsistencyMethodProj), INTENT(INOUT) :: consistency

		INTEGER :: i
		DOUBLE PRECISION :: b

		DOUBLE PRECISION :: eHF
		DOUBLE PRECISION :: kinetic_energy, kinetic_CM_energy
		DOUBLE PRECISION :: local_energy_BB, exchange_energy_BB ! Brink-Booker
		DOUBLE PRECISION :: local_energy_Coulomb, exchange_energy_Coulomb
		DOUBLE PRECISION :: dd_energy ! Density Dependent
		DOUBLE PRECISION :: ls_energy

		DOUBLE PRECISION :: pairing
		DOUBLE PRECISION :: pairing_Ek
		DOUBLE PRECISION :: pairing_BB ! BrinkBooker
		DOUBLE PRECISION :: pairing_Coulomb
		DOUBLE PRECISION :: pairing_LS !

		DOUBLE PRECISION :: N, Z

		b = Nucleus_get_b(consistency%density%nucleus)

		WRITE(*,'(/,5X,"SUMMARY OF THE RUN - ENERGIES")')
		WRITE(*,'(5X,"=============================")')

		PRINT *
		PRINT "(A50,F20.10)", "Longitud del oscilador:", b
		PRINT "(70A1)", ("-", i = 1, 70)

		kinetic_energy = SelfConsistencyMethodProj_get_Ek(consistency)
		kinetic_CM_energy = SelfConsistencyMethodProj_get_EkCM(consistency)
		local_energy_BB = SelfConsistencyMethodProj_get_LocalBBEnergy(consistency)
 		local_energy_Coulomb = SelfConsistencyMethodProj_get_LocalCoulombEnergy(consistency)
		exchange_energy_BB = SelfConsistencyMethodProj_get_ExchangeBBEnergy(consistency)
		exchange_energy_Coulomb = SelfConsistencyMethodProj_get_ExchangeCoulombEnergy(consistency)
		dd_energy = SelfConsistencyMethodProj_get_DDEnergy(consistency)
		ls_energy = SelfConsistencyMethodProj_get_LSEnergy(consistency)

		PRINT "(A50,F20.10)", "Energia cinetica:", kinetic_energy
		PRINT "(A50,F20.10)", "Energia cinetica del centro de masas:", kinetic_CM_energy
		PRINT "(A50,F20.10)", "Energia local de Brink-Booker:", local_energy_BB
		PRINT "(A50,F20.10)", "Energia local de Coulomb:", local_energy_Coulomb
		PRINT "(A50,F20.10)", "Energia de intercambio de Brink-Booker:", exchange_energy_BB
		PRINT "(A50,F20.10)", "Energia de intercambio de Coulomb:", exchange_energy_Coulomb
		PRINT "(A50,F20.10)", "Energia dependiente de la densidad:", dd_energy
		PRINT "(A50,F20.10)", "Energia de spin-orbita:", ls_energy
		PRINT "(70A1)", ("-", i = 1, 70)

		eHF = kinetic_energy + kinetic_CM_energy &
			+ local_energy_BB + exchange_energy_BB &
			+ local_energy_Coulomb + exchange_energy_Coulomb &
			+ dd_energy + ls_energy

		PRINT "(A50,F20.10)", "Energia total de Brink-Booker:", local_energy_BB + exchange_energy_BB
		PRINT "(A50,F20.10)", "Energia total de Coulomb:", local_energy_Coulomb + exchange_energy_Coulomb
		PRINT "(A50,F20.10)", "Energia total (Hartree-Fock):", eHF
		PRINT "(70A1)", ("-", i = 1, 70)

		pairing_Ek = SelfConsistencyMethodProj_get_EkPairing(consistency)
		pairing_BB = SelfConsistencyMethodProj_get_BBPairing(consistency)
		pairing_Coulomb = SelfConsistencyMethodProj_get_CoulombPairing(consistency)
		pairing_LS = SelfConsistencyMethodProj_get_LSPairing(consistency)

		PRINT "(A50,E20.10)", "Apareamiento del centro de masas:", pairing_Ek
		PRINT "(A50,E20.10)", "Apareamiento de Brink-Booker:", pairing_BB
		PRINT "(A50,E20.10)", "Apareamiento de Coulomb:", pairing_Coulomb
		PRINT "(A50,E20.10)", "Apareamiento de spin-orbita:", pairing_LS
		PRINT "(70A1)", ("-", i = 1, 70)

		pairing = pairing_BB + pairing_Coulomb + pairing_LS + pairing_Ek

		PRINT "(A50,E20.10)", "Apareamiento total:", pairing
		PRINT "(70A1)", ("-", i = 1, 70)

		PRINT "(A50,F20.10)", "TOTAL:", eHF + pairing

		PRINT *
		PRINT "(A40,A15,A15)", "Neutrones", "Protones", "Total"

		N = consistency%density%nucleus%actual_np(1)
		Z = consistency%density%nucleus%actual_np(0)

		PRINT "(A25,F15.5,F15.5,F15.5)", "Particulas:", N, Z, N + Z
		PRINT "(A25,F15.5,F15.5)", "Potenciales quimicos:", &
			consistency%density%nucleus%lambda_np(1), &
			consistency%density%nucleus%lambda_np(0)
		PRINT "(A25,F15.5,F15.5,F15.5)", "Radio:", &
			SQRT(consistency%density%nucleus%actual_R2(1)), &
			SQRT(consistency%density%nucleus%actual_R2(0)), &
			SQRT(Nucleus_get_actual_R2(consistency%density%nucleus))
		PRINT "(A25,F15.5,F15.5,F15.5)", "Multiplicadores y error:", &
			consistency%density%nucleus%lambda_R2(1), &
			consistency%density%nucleus%lambda_R2(0), &
			Nucleus_get_actual_R2(consistency%density%nucleus) - consistency%density%nucleus%R2
! D.showSpatialDistribution(NEU)
!TODO		CALL Nucleus_show_ExperimentalData(consistency%density%nucleus)
		RETURN
	END SUBROUTINE SelfConsistencyMethodProj_show_Status

	SUBROUTINE SelfConsistencyMethodProj_del(consistency)
		TYPE (SelfConsistencyMethodProj), INTENT(INOUT) :: consistency

		CALL SymKineticEnergy2Body_del(consistency%vEkCMph, consistency%vEkCMpp)

		CALL SymVBBph_del(consistency%vBBph)
		CALL SymVBBpp_del(consistency%vBBpp)
		CALL SymVCph_del(consistency%vCph)
		CALL SymVCpp_del(consistency%vCpp)
		CALL SymVLSph_del(consistency%vLSph)
		CALL SymVLSpp_del(consistency%vLSpp)
		CALL SymGDDph_del(consistency%gDDph)

		NULLIFY(consistency%density)

		RETURN
	END SUBROUTINE SelfConsistencyMethodProj_del

END MODULE selfc_proj
