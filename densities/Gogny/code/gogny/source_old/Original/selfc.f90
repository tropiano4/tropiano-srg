MODULE selfc

	USE input
	USE nucleus
	USE symden
	USE symke2b
	USE symvbb
	USE symvc
	USE symvls
	USE symgdd

	IMPLICIT NONE

	TYPE SelfConsistencyMethod
		TYPE (SymKineticEnergy2Body) vEkCMph
		TYPE (SymVBBph) vBBph
		TYPE (SymVBBpp) vBBpp
		TYPE (SymVCph) vCph
		TYPE (SymVCpp) vCpp
		TYPE (SymVLSph) vLSph
		TYPE (SymVLSpp) vLSpp
		TYPE (SymGDDph) gDDph
		TYPE (SymDensity), POINTER :: density
	END TYPE

CONTAINS

	SUBROUTINE SelfConsistencyMethod_new(consistency, density)
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency
		TYPE (SymDensity), TARGET, INTENT(IN) :: density

		DOUBLE PRECISION b

		b = Nucleus_get_b(density%nucleus)
		CALL SymKineticEnergy2Body_new(consistency%vEkCMph)
		CALL SymVBBph_new(consistency%vBBph, b)
		CALL SymVBBpp_new(consistency%vBBpp, b)
		CALL SymVCph_new(consistency%vCph)
		CALL SymVCpp_new(consistency%vCpp)
		CALL SymVLSph_new(consistency%vLSph)
		CALL SymVLSpp_new(consistency%vLSpp)
		CALL SymGDDph_new(consistency%gDDph)
		consistency%density => density
		RETURN
	END SUBROUTINE SelfConsistencyMethod_new

	SUBROUTINE SelfConsistencyMethod_store_eHFB(consistency)
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency

		DOUBLE PRECISION eHF
		DOUBLE PRECISION kinetic_energy, kinetic_CM_energy
		DOUBLE PRECISION local_energy_BB, exchange_energy_BB ! Brink-Booker
		DOUBLE PRECISION local_energy_Coulomb, exchange_energy_Coulomb
		DOUBLE PRECISION dd_energy ! Density Dependent
		DOUBLE PRECISION ls_energy

		DOUBLE PRECISION pairing
		DOUBLE PRECISION pairing_BB ! BrinkBooker
		DOUBLE PRECISION pairing_Coulomb
		DOUBLE PRECISION pairing_LS !

		kinetic_energy          = SelfConsistencyMethod_get_Ek(consistency)
		kinetic_CM_energy       = SelfConsistencyMethod_get_EkCM(consistency)
		local_energy_BB         = SelfConsistencyMethod_get_LocalBBEnergy(consistency)
 		local_energy_Coulomb    = SelfConsistencyMethod_get_LocalCoulombEnergy(consistency)
		exchange_energy_BB      = SelfConsistencyMethod_get_ExchangeBBEnergy(consistency)
		exchange_energy_Coulomb = SelfConsistencyMethod_get_ExchangeCoulombEnergy(consistency)
		dd_energy               = SelfConsistencyMethod_get_DDEnergy(consistency)
		ls_energy               = SelfConsistencyMethod_get_LSEnergy(consistency)

		eHF = kinetic_energy + kinetic_CM_energy &
			+ local_energy_BB + exchange_energy_BB &
			+ local_energy_Coulomb + exchange_energy_Coulomb &
			+ dd_energy + ls_energy

		pairing_BB      = SelfConsistencyMethod_get_BBPairing(consistency)
		pairing_Coulomb = SelfConsistencyMethod_get_CoulombPairing(consistency)
		pairing_LS      = SelfConsistencyMethod_get_LSPairing(consistency)

		pairing = pairing_BB + pairing_Coulomb + pairing_LS

!	print *, "eHF & pairing =", eHF, pairing
!	print *, "kinetic_energy=", kinetic_energy
!	print *, "kinetic_CM_energy=", kinetic_CM_energy
!	print *, "local_energy_BB=", local_energy_BB
!	print *, "local_energy_Coulomb=", local_energy_Coulomb
!	print *, "exchange_energy_BB=", exchange_energy_BB
!	print *, "exchange_energy_Coulomb=", exchange_energy_Coulomb
!	print *, "dd_energy=", dd_energy
!	print *, "ls_energy=", ls_energy
!	print *, "pairing_BB=", pairing_BB
!	print *, "pairing_Coulomb=", pairing_Coulomb
!	print *, "pairing_LS=", pairing_LS
		consistency%density%nucleus%eHFB = eHF + pairing
		consistency%density%nucleus%pairing = pairing
		RETURN
	END SUBROUTINE SelfConsistencyMethod_store_eHFB

	! Energía cinética
	FUNCTION SelfConsistencyMethod_get_Ek(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_Ek
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		INTEGER A
		DOUBLE PRECISION b
		DOUBLE PRECISION factor
		INTEGER ta
		DOUBLE PRECISION, DIMENSION(0:1) :: Ek

		A = Nucleus_get_A(consistency%density%nucleus)
		b = Nucleus_get_b(consistency%density%nucleus)
		factor = (1.0 - (1.0 / A)) / (b ** 2.0)
		DO ta = 0, 1
			Ek(ta) = factor * (EkField * consistency%density%field%rho%p(ta))
		END DO
		SelfConsistencyMethod_get_Ek = Ek(0) + Ek(1)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_Ek

	FUNCTION SelfConsistencyMethod_get_EkCM(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_EkCM
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		INTEGER A
		DOUBLE PRECISION b
		TYPE (SymHartreeFockField) HF_Gamma

		CALL SymHartreeFockField_new(HF_Gamma)
		CALL SymKineticEnergy2Body_get_Gamma(HF_Gamma, consistency%vEkCMph, consistency%density%field%rho)

		A = Nucleus_get_A(consistency%density%nucleus)
		b = Nucleus_get_b(consistency%density%nucleus)
		IF ((A .LE. 1) .OR. (b .LE. 0.1)) STOP "Abortado"
		SelfConsistencyMethod_get_EkCM = (0.5 / (A * (b ** 2))) * (consistency%density%field%rho * HF_Gamma)

		CALL SymHartreeFockField_del(HF_Gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_EkCM

	! Energía local de Brink-Booker
	FUNCTION SelfConsistencyMethod_get_LocalBBEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_LocalBBEnergy
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) local_gamma

		CALL SymHartreeFockField_new(local_gamma)
		CALL SymVBBph_get_LocalGamma(local_gamma, consistency%vBBph, consistency%density%field%rho)

		SelfConsistencyMethod_get_LocalBBEnergy = 0.5 * (consistency%density%field%rho * local_gamma)

		CALL SymHartreeFockField_del(local_gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_LocalBBEnergy

	! Energía de intercambio de Brink-Booker
	FUNCTION SelfConsistencyMethod_get_ExchangeBBEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_ExchangeBBEnergy
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) exchange_gamma

		CALL SymHartreeFockField_new(exchange_gamma)
		CALL SymVBBph_get_ExchangeGamma(exchange_gamma, consistency%vBBph, consistency%density%field%rho)

		SelfConsistencyMethod_get_ExchangeBBEnergy = 0.5 * (consistency%density%field%rho * exchange_gamma)

		CALL SymHartreeFockField_del(exchange_gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_ExchangeBBEnergy

	! Energía local de Coulomb
	FUNCTION SelfConsistencyMethod_get_LocalCoulombEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_LocalCoulombEnergy
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		DOUBLE PRECISION b
		TYPE (SymHartreeFockField) local_gamma

		CALL SymHartreeFockField_new(local_gamma)
		CALL SymVCph_get_LocalGamma(local_gamma, consistency%vCph, consistency%density%field%rho)

		b = Nucleus_get_b(consistency%density%nucleus)
		SelfConsistencyMethod_get_LocalCoulombEnergy = (0.5 / b) * (consistency%density%field%rho * local_gamma)

		CALL SymHartreeFockField_del(local_gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_LocalCoulombEnergy

	! Energía de intercambio de Coulomb
	FUNCTION SelfConsistencyMethod_get_ExchangeCoulombEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_ExchangeCoulombEnergy
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		DOUBLE PRECISION b
		TYPE (SymHartreeFockField) exchange_gamma

		CALL SymHartreeFockField_new(exchange_gamma)
		CALL SymVCph_get_ExchangeGamma(exchange_gamma, consistency%vCph, consistency%density%field%rho)

		b = Nucleus_get_b(consistency%density%nucleus)
		SelfConsistencyMethod_get_ExchangeCoulombEnergy = (0.5 / b) * (consistency%density%field%rho * exchange_gamma)

		CALL SymHartreeFockField_del(exchange_gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_ExchangeCoulombEnergy

	! Energía dependiente de la densidad
	FUNCTION SelfConsistencyMethod_get_DDEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_DDEnergy
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency

		DOUBLE PRECISION b

		b = Nucleus_get_b(consistency%density%nucleus)
		CALL SymGDDph_make_DD(consistency%gDDph, consistency%density%field%rho)
		SelfConsistencyMethod_get_DDEnergy = SymGDDph_get_edd(consistency%gDDph) / (b ** 4.0)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_DDEnergy

	FUNCTION SelfConsistencyMethod_get_LSEnergy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_LSEnergy
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		DOUBLE PRECISION b
		TYPE (SymHartreeFockField) HF_Gamma

		CALL SymHartreeFockField_new(HF_Gamma)
		CALL SymVLSph_get_Gamma(HF_Gamma, consistency%vLSph, consistency%density%field%rho)

		b = Nucleus_get_b(consistency%density%nucleus)
		SelfConsistencyMethod_get_LSEnergy = (0.5 / (b ** 5.0)) * (consistency%density%field%rho * HF_Gamma)

		CALL SymHartreeFockField_del(HF_Gamma)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_LSEnergy

	! Apareamiento de Brink-Booker
	FUNCTION SelfConsistencyMethod_get_BBPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_BBPairing
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		TYPE (SymHartreeFockField) vbb_field

		CALL SymHartreeFockField_new(vbb_field)
		CALL SymVBBpp_product(vbb_field, consistency%vBBpp, consistency%density%field%kap)

		SelfConsistencyMethod_get_BBPairing = 0.5 * (consistency%density%field%kap * vbb_field)

		CALL SymHartreeFockField_del(vbb_field)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_BBPairing

	! Apareamiento de Coulomb
	FUNCTION SelfConsistencyMethod_get_CoulombPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_CoulombPairing
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		DOUBLE PRECISION b
		TYPE (SymHartreeFockField) delta

		CALL SymHartreeFockField_new(delta)
		CALL SymVCpp_get_Delta(delta, consistency%vCpp, consistency%density%field%kap)

		b = Nucleus_get_b(consistency%density%nucleus)
		SelfConsistencyMethod_get_CoulombPairing = (0.5 / b) * (consistency%density%field%kap * delta)

		CALL SymHartreeFockField_del(delta)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_CoulombPairing

	FUNCTION SelfConsistencyMethod_get_LSPairing(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_get_LSPairing
		TYPE (SelfConsistencyMethod), INTENT(IN) :: consistency

		DOUBLE PRECISION b
		TYPE (SymHartreeFockField) delta

		CALL SymHartreeFockField_new(delta)
		CALL SymVLSpp_get_Delta(delta, consistency%vLSpp, consistency%density%field%kap)

		b = Nucleus_get_b(consistency%density%nucleus)
		SelfConsistencyMethod_get_LSPairing = (0.5 / (b ** 5.0)) * (consistency%density%field%kap * delta)

		CALL SymHartreeFockField_del(delta)
		RETURN
	END FUNCTION SelfConsistencyMethod_get_LSPairing

	FUNCTION SelfConsistencyMethod_accuracy(consistency)
		DOUBLE PRECISION SelfConsistencyMethod_accuracy
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency

		DOUBLE PRECISION old_eHFB

		old_eHFB = consistency%density%nucleus%eHFB
		CALL SelfConsistencyMethod_store_eHFB(consistency)
		SelfConsistencyMethod_accuracy = ABS(old_eHFB - consistency%density%nucleus%eHFB)
		RETURN
	END FUNCTION SelfConsistencyMethod_accuracy

	SUBROUTINE SelfConsistencyMethod_show_Status(consistency)
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency

		INTEGER i
		DOUBLE PRECISION b

		DOUBLE PRECISION eHF
		DOUBLE PRECISION kinetic_energy, kinetic_CM_energy
		DOUBLE PRECISION local_energy_BB, exchange_energy_BB ! Brink-Booker
		DOUBLE PRECISION local_energy_Coulomb, exchange_energy_Coulomb
		DOUBLE PRECISION dd_energy ! Density Dependent
		DOUBLE PRECISION ls_energy

		DOUBLE PRECISION pairing
		DOUBLE PRECISION pairing_BB ! BrinkBooker
		DOUBLE PRECISION pairing_Coulomb
		DOUBLE PRECISION pairing_LS !

		DOUBLE PRECISION N, Z

		b = Nucleus_get_b(consistency%density%nucleus)
		PRINT *
		PRINT "(A50,F20.10)", "Longitud del oscilador:", b
		PRINT "(70A1)", ("-", i = 1, 70)

		kinetic_energy = SelfConsistencyMethod_get_Ek(consistency)
		kinetic_CM_energy = SelfConsistencyMethod_get_EkCM(consistency)
		local_energy_BB = SelfConsistencyMethod_get_LocalBBEnergy(consistency)
 		local_energy_Coulomb = SelfConsistencyMethod_get_LocalCoulombEnergy(consistency)
		exchange_energy_BB = SelfConsistencyMethod_get_ExchangeBBEnergy(consistency)
		exchange_energy_Coulomb = SelfConsistencyMethod_get_ExchangeCoulombEnergy(consistency)
		dd_energy = SelfConsistencyMethod_get_DDEnergy(consistency)
		ls_energy = SelfConsistencyMethod_get_LSEnergy(consistency)
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
		PRINT "(A50,F20.10)", "Enertia total de Brink-Booker:", local_energy_BB + exchange_energy_BB
		PRINT "(A50,F20.10)", "Enertia total de Coulomb:", local_energy_Coulomb + exchange_energy_Coulomb
		PRINT "(A50,F20.10)", "Energia total (Hartree-Fock):", eHF
		PRINT "(70A1)", ("-", i = 1, 70)

		pairing_BB = SelfConsistencyMethod_get_BBPairing(consistency)
		pairing_Coulomb = SelfConsistencyMethod_get_CoulombPairing(consistency)
		pairing_LS = SelfConsistencyMethod_get_LSPairing(consistency)
		PRINT "(A50,E20.10)", "Apareamiento de Brink-Booker:", pairing_BB
		PRINT "(A50,E20.10)", "Apareamiento de Coulomb:", pairing_Coulomb
		PRINT "(A50,E20.10)", "Apareamiento de spin-orbita:", pairing_LS
		PRINT "(70A1)", ("-", i = 1, 70)

		pairing = pairing_BB + pairing_Coulomb + pairing_LS
		PRINT "(A50,E20.10)", "Apareamiento total:", pairing
		PRINT "(70A1)", ("-", i = 1, 70)

		PRINT "(A50,F20.10)", "TOTAL:", eHF + pairing

		PRINT *
		PRINT "(A40,A15,A15)", "Neutrones", "Protones", "Total"
		N = consistency%density%nucleus%actual_np(0)
		Z = consistency%density%nucleus%actual_np(1)
		PRINT "(A25,F15.5,F15.5,F15.5)", "Particulas:", N, Z, N + Z
		PRINT "(A25,F15.5,F15.5)", "Potenciales quimicos:", &
			consistency%density%nucleus%lambda_np(0), &
			consistency%density%nucleus%lambda_np(1)
		PRINT "(A25,F15.5,F15.5,F15.5)", "Radio:", &
			SQRT(consistency%density%nucleus%actual_R2(0)), &
			SQRT(consistency%density%nucleus%actual_R2(1)), &
			SQRT(Nucleus_get_actual_R2(consistency%density%nucleus))
		PRINT "(A25,F15.5,F15.5,F15.5)", "Multiplicadores y error:", &
			consistency%density%nucleus%lambda_R2(0), &
			consistency%density%nucleus%lambda_R2(1), &
			Nucleus_get_actual_R2(consistency%density%nucleus) - consistency%density%nucleus%R2
! D.showSpatialDistribution(NEU)
!TODO		CALL Nucleus_show_ExperimentalData(consistency%density%nucleus)
		RETURN
	END SUBROUTINE SelfConsistencyMethod_show_Status

	SUBROUTINE SelfConsistencyMethod_del(consistency)
		TYPE (SelfConsistencyMethod), INTENT(INOUT) :: consistency

		CALL SymKineticEnergy2Body_del(consistency%vEkCMph)
		CALL SymVBBph_del(consistency%vBBph)
		CALL SymVBBpp_del(consistency%vBBpp)
		CALL SymVCph_del(consistency%vCph)
		CALL SymVCpp_del(consistency%vCpp)
		CALL SymVLSph_del(consistency%vLSph)
		CALL SymVLSpp_del(consistency%vLSpp)
		CALL SymGDDph_del(consistency%gDDph)
		NULLIFY(consistency%density)
		RETURN
	END SUBROUTINE SelfConsistencyMethod_del

END MODULE selfc
