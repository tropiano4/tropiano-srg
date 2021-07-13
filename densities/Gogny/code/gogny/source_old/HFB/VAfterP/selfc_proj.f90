 MODULE selfc_proj

	USE input
	USE nucleus
	USE symden_proj
	USE symke2b
	USE symvbb
	USE symvc
	USE symvls
	USE symgdd_proj
	USE energy_proj

	IMPLICIT NONE

	TYPE SelfConsistencyMethodProj
		TYPE (SymKineticEnergy2Body) vEkCMph
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

		pairing_BB      = SelfConsistencyMethodProj_get_BBPairing(consistency)
		pairing_Coulomb = SelfConsistencyMethodProj_get_CoulombPairing(consistency)
		pairing_LS      = SelfConsistencyMethodProj_get_LSPairing(consistency)

		pairing = pairing_BB + pairing_Coulomb + pairing_LS

		consistency%density%nucleus%eHFB = eHF + pairing
		consistency%density%nucleus%pairing = pairing
		
		RETURN
	END SUBROUTINE SelfConsistencyMethodProj_store_eHFB

	! One-body Kinetic Energy
	FUNCTION SelfConsistencyMethodProj_get_Ek(consistency)
		DOUBLE PRECISION SelfConsistencyMethodProj_get_Ek
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency

		INTEGER A
		DOUBLE PRECISION :: factor, b
		INTEGER ta
		DOUBLE PRECISION, DIMENSION(0:1) :: Ek

		A = Nucleus_get_A(consistency%density%nucleus)
                b = Nucleus_get_b(consistency%density%nucleus)
		
                factor = (1.0 - (1.0 / A))
		
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
		
		SelfConsistencyMethodProj_get_EkCM = (0.5 / A) * (HF_Gamma * consistency%density%field%rho)

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
		
		SelfConsistencyMethodProj_get_LocalBBEnergy = 0.5 * (HF_Gamma * consistency%density%field%rho)

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
		
		SelfConsistencyMethodProj_get_ExchangeBBEnergy = 0.5 * (HF_Gamma * consistency%density%field%rho)

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

		SelfConsistencyMethodProj_get_LocalCoulombEnergy = 0.5 * (HF_Gamma * consistency%density%field%rho)

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

		SelfConsistencyMethodProj_get_ExchangeCoulombEnergy = 0.5 * (HF_Gamma * consistency%density%field%rho)

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
		
		SelfConsistencyMethodProj_get_DDEnergy = SymGDDph_get_edd(consistency%gDDph)
		
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

		SelfConsistencyMethodProj_get_LSEnergy = 0.5* (HF_Gamma * consistency%density%field%rho)

		CALL SymHartreeFockFieldProj_del(HF_Gamma)
		
		RETURN
	END FUNCTION SelfConsistencyMethodProj_get_LSEnergy

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

		IF (HFOnly .EQ. 0) THEN
			SelfConsistencyMethodProj_get_BBPairing = 0.5 * (consistency%density%field%kap * delta)
		ELSE
			SelfConsistencyMethodProj_get_BBPairing = 0.0
		END IF

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

		IF (HFOnly .EQ. 0) THEN
			SelfConsistencyMethodProj_get_CoulombPairing = 0.5 * (consistency%density%field%kap * delta)
		ELSE
			SelfConsistencyMethodProj_get_CoulombPairing = 0.0
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

		IF (HFOnly .EQ. 0) THEN
			SelfConsistencyMethodProj_get_LSPairing = 0.5 * (consistency%density%field%kap * delta)
		ELSE
			SelfConsistencyMethodProj_get_LSPairing = 0.0
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
			
		PRINT "(A50,F20.10)", "Enertia total de Brink-Booker:", local_energy_BB + exchange_energy_BB
		PRINT "(A50,F20.10)", "Enertia total de Coulomb:", local_energy_Coulomb + exchange_energy_Coulomb
		PRINT "(A50,F20.10)", "Energia total (Hartree-Fock):", eHF
		PRINT "(70A1)", ("-", i = 1, 70)

		pairing_BB = SelfConsistencyMethodProj_get_BBPairing(consistency)
		pairing_Coulomb = SelfConsistencyMethodProj_get_CoulombPairing(consistency)
		pairing_LS = SelfConsistencyMethodProj_get_LSPairing(consistency)
		
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
	END SUBROUTINE SelfConsistencyMethodProj_del
	
	
	
	
	
	
	
	!-------------------------------------------------------------------------------!
	!										!
	!   	CALCULATING THE PROJECTED HFB MATRIX: MEAN-FIELD (Gamma)		!
	!										!
	!  Calculate Gamma_{ta,tb}(phi)							!
	!										!
	!  The field is obtained with the matrix elements corresponding to isospin ta	!
	!  and the density corresponding to isospin tb. 				!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGammaTensorVAP_get_GammaPhi(HF_out, genden_gauge, genden_proj, isoCoupling, AngleGauge, consistency)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SelfConsistencyMethodProj), INTENT(OUT) :: consistency
		TYPE (SymGenDensityGaugeProj), INTENT(IN) :: genden_gauge
		
		DOUBLE PRECISION, INTENT(IN) :: AngleGauge
		INTEGER, INTENT(IN) :: isoCoupling
		
		TYPE (SymHartreeFockFieldProj) :: ekcm_field, vbb_field, vc_field, vls_field, gdd_field
		TYPE (SymHartreeFockFieldProj) :: field1, field2, rho_field

		TYPE (SymGenDensityProj) :: genden_proj
		TYPE (SymDensityProj) :: density, density_proj		

		COMPLEX :: factor, Trace_VAP
		DOUBLE PRECISION :: A
		INTEGER :: tb, ta
		
		A = Nucleus_get_A(consistency%density%nucleus)
		
		CALL SymHartreeFockFieldProj_new(ekcm_field)
		CALL SymHartreeFockFieldProj_new(vbb_field)
		CALL SymHartreeFockFieldProj_new(vc_field)
		CALL SymHartreeFockFieldProj_new(vls_field)
		CALL SymHartreeFockFieldProj_new(gdd_field)
		
		CALL SymHartreeFockFieldProj_new(field1)
		CALL SymHartreeFockFieldProj_new(field2)
		CALL SymHartreeFockFieldProj_new(rho_field)
		
		! Create new gauge-dependent density "density": rho(phi)
		
		!CALL SymDensity_new_GenDensityProj10(density, genden_gauge)
		!CALL SymDensity_new_GenDensityProj(density_proj, genden_proj)
		
		!CALL SymGenDensityHFProj_assign_VAP(rho_field, genden_gauge%rho, 0)
		CALL SymGenDensityHFProj_assign2(rho_field, genden_gauge%rho)
		
		DO ta = 0, 1
		
						tb = 1-ta
			IF (isoCoupling .EQ. 0) tb = ta
			
			!CALL SymKineticEnergy2Body_product(field1, consistency%vEkCMph, density%field%rho, ta, tb, 1)
			CALL SymKineticEnergy2Body_product(field1, consistency%vEkCMph, rho_field, ta, tb, 1)
			CALL SymHartreeFockFieldProj_product(ekcm_field, CMPLX(1.0/A, 0), field1)
		
			! Mean-field - Brink-Boeker term
		
			!CALL SymVBBph_get_LocalGamma(field1, consistency%vBBph, density%field%rho, ta, tb, 1)
			CALL SymVBBph_get_LocalGamma(field1, consistency%vBBph, rho_field, ta, tb, 1)
			!CALL SymVBBph_get_ExchangeGamma(field2, consistency%vBBph, density%field%rho, ta, tb, 1)
			CALL SymVBBph_get_ExchangeGamma(field2, consistency%vBBph, rho_field, ta, tb, 1)
		
			CALL SymHartreeFockFieldProjIso_add(vbb_field, field2, field1, ta)

			! Mean-field - Coulomb potential
			!CALL SymVCph_get_Gamma(vc_field, consistency%vCph, density%field%rho, ta, tb, 1)
			CALL SymVCph_get_Gamma(vc_field, consistency%vCph, rho_field, ta, tb, 1)

			! Mean-field - Spin-orbit term
			!CALL SymVLSph_get_Gamma(vls_field, consistency%vLSph, density%field%rho, ta, tb, 1)
			CALL SymVLSph_get_Gamma(vls_field, consistency%vLSph, rho_field, ta, tb, 1)

			! Mean-field - Density-dependent term. This term needs to be calculated with the projected density, contrary
			!              all the others.
			!IF (tb .EQ. ta) THEN
			!	CALL SymGDDph_update(gdd_field, consistency%gDDph, density_proj%field%rho, AngleGauge)
			!END IF

			! Total Mean-field = Sum of all the preceding terms

			CALL SymHartreeFockFieldProjIso_add(field1, vls_field, vbb_field, ta)
			CALL SymHartreeFockFieldProjIso_add(field2, vc_field, field1, ta)
			!CALL SymHartreeFockFieldProjIso_add(field1, gdd_field, field2, ta)
			CALL SymHartreeFockFieldProjIso_add(HF_out, ekcm_field, field2, ta)
				
		END DO
		
		CALL SymHartreeFockFieldProj_assign(HF_out, vbb_field)
		
		CALL SymHartreeFockFieldProj_del(ekcm_field)
		CALL SymHartreeFockFieldProj_del(vbb_field)
		CALL SymHartreeFockFieldProj_del(vc_field)
		CALL SymHartreeFockFieldProj_del(vls_field)
		CALL SymHartreeFockFieldProj_del(gdd_field)

		CALL SymHartreeFockFieldProj_del(rho_field)
		CALL SymHartreeFockFieldProj_del(field1)
		CALL SymHartreeFockFieldProj_del(field2)

		RETURN
	END SUBROUTINE SymGammaTensorVAP_get_GammaPhi
	
	!-------------------------------------------------------------------------------!
	!										!
	!   	CALCULATING THE PROJECTED HFB MATRIX: MEAN-FIELD (Gamma)		!
	!										!
	!  Calculate Gamma_{DD, ta,tb}(phi)						!
	!										!
	!  The field is obtained with the matrix elements corresponding to isospin ta	!
	!  and the density corresponding to isospin tb. 				!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGammaTensorVAP_get_GammaDDPhi(HF_out, Mfactor, CoeffsXY, genden_gauge, C_Matrix, Deriv_Y, NGauge, IndexGauge, ta)
		TYPE (SymHartreeFockFieldProj), INTENT(OUT) :: HF_out
		TYPE (ProjectionCoeffs), INTENT(IN) :: CoeffsXY
		TYPE (SymHartreeFockFieldProj), DIMENSION(:), POINTER :: Mfactor
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix, Deriv_Y
		
		INTEGER, INTENT(IN) :: ta, IndexGauge, NGauge
		
		TYPE (SymHartreeFockFieldProj) :: field1, HF_1, HF_2, dY
		TYPE (SymDensityProj), DIMENSION(:), POINTER :: density_gauge
		TYPE (SymDensityProj) :: density, density_proj
		
		COMPLEX :: Trace_1
		DOUBLE PRECISION :: pi, AngleGauge 
		INTEGER :: IndexGauge_tau_prime
			
		CALL SymHartreeFockFieldProj_new(field1)
		CALL SymHartreeFockFieldProj_new(HF_1)
		CALL SymHartreeFockFieldProj_new(HF_2)
		CALL SymHartreeFockFieldProj_new(dY)
		
		pi = 4.0*ATAN(1.0)
		
		! Transform the density into the proper form before to do anything 
		! (required because when "converting" densities, both isospins are treated
		! simultaneously)
		
		ALLOCATE(density_gauge(1:NGauge))
		
		DO IndexGauge_tau_prime = 1, NGauge
			CALL SymDensity_new_GenDensityProj_only_rho(density_gauge(IndexGauge_tau_prime), genden_gauge(IndexGauge_tau_prime))			
		END DO
				
		DO IndexGauge_tau_prime = 1, NGauge
		
			AngleGauge = pi*REAL(IndexGauge_tau_prime)/REAL(NGauge)

			! Create new gauge-dependent densities: rho(phi) (density_gauge) and dy/d_rho (dY)
			
			CALL SymGenDensityHFProj_assign_VAP(dY, Deriv_Y(IndexGauge_tau_prime), 0)
							
			Trace_1 = SymHartreeFockFieldProj_product_iso(Mfactor(IndexGauge), density_gauge(IndexGauge_tau_prime)%field%rho, ta)
				
			CALL SymHartreeFockFieldProj_scalar_iso(field1, 0.5*Trace_1, dY, ta)
				
			CALL SymGammaTensorVAP_get_GammaDDProj(genden_gauge, C_Matrix, ta, AngleGauge, IndexGauge_tau_prime, &
													Mfactor(IndexGauge), HF_1)
											
			CALL SymHartreeFockFieldProj_scalar_iso(HF_2, 0.5*CoeffsXY%y_l(IndexGauge_tau_prime, ta), HF_1, ta)
					
			CALL SymHartreeFockFieldProjIso_add(HF_1, field1, HF_2, ta)
			CALL SymHartreeFockFieldProjIso_add(HF_out, HF_out, HF_1, ta)
					
		END DO

		CALL SymHartreeFockFieldProj_del(field1)
		CALL SymHartreeFockFieldProj_del(HF_1)
		CALL SymHartreeFockFieldProj_del(HF_2)
		CALL SymHartreeFockFieldProj_del(dY)
		
		RETURN
	END SUBROUTINE SymGammaTensorVAP_get_GammaDDPhi
	
	!-------------------------------------------------------------------------------!
	!										!
	!   	CALCULATING THE PROJECTED HFB MATRIX: MEAN-FIELD (Gamma)		!
	!										!
	!  Calculate Gamma_{ta,tb}(phi)							!
	!										!
	!  The field is obtained with the matrix elements corresponding to isospin ta	!
	!  and the density corresponding to isospin tb. 				!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGammaTensorVAP_get_DeltaPhi(HF_out, genden_gauge, consistency, TypeKap)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SelfConsistencyMethodProj), INTENT(IN) :: consistency
		TYPE (SymGenDensityGaugeProj), INTENT(IN) :: genden_gauge
		
		INTEGER, INTENT(IN) :: TypeKap
		
		TYPE (SymHartreeFockFieldProj) :: kap_field
		TYPE (SymHartreeFockFieldProj) :: vbb_field, vc_field, vls_field, Delta
		TYPE (SymDensityProj) :: density
				
		INTEGER :: ta

		CALL SymHartreeFockFieldProj_new(vbb_field)
		CALL SymHartreeFockFieldProj_new(vc_field)
		CALL SymHartreeFockFieldProj_new(vls_field)
		CALL SymHartreeFockFieldProj_new(Delta)
		
		CALL SymHartreeFockFieldProj_new(kap_field)
		
		! Create new gauge-dependent densities kappa01 and kappa10
		IF (TypeKap .EQ. 0) THEN
			!CALL SymDensity_new_GenDensityProj01(density, genden_gauge)
			CALL SymGenDensityHFProj_assign_VAP(kap_field, genden_gauge%kap01, 0)
		ELSE
			!CALL SymDensity_new_GenDensityProj10(density, genden_gauge)
			CALL SymGenDensityHFProj_assign_VAP(kap_field, genden_gauge%kap10, 0)
		END IF
		
		DO ta = 0, 1
				
			! Pairing - Brink-Boker term
			!CALL SymVBBpp_get_Delta(vbb_field, consistency%vBBpp, density%field%kap, ta, ta)
			CALL SymVBBpp_get_Delta(vbb_field, consistency%vBBpp, kap_field, ta, ta)

			! Pairing - Coulomb term
			!CALL SymVCpp_get_Delta(vc_field, consistency%vCpp, density%field%kap, ta, ta)
			CALL SymVCpp_get_Delta(vc_field, consistency%vCpp, kap_field, ta, ta)

			! Pairing - Spin-orbit term
			!CALL SymVLSpp_get_Delta(vls_field, consistency%vLSpp, density%field%kap, ta, ta)
			CALL SymVLSpp_get_Delta(vls_field, consistency%vLSpp, kap_field, ta, ta)

		END DO
		
		! Total Pairing = Sum of all the preceding terms

		CALL SymHartreeFockFieldProj_add(Delta, vc_field, vls_field)
		CALL SymHartreeFockFieldProj_add(HF_out, vbb_field, Delta)
			
		CALL SymHartreeFockFieldProj_del(kap_field)
		
		CALL SymHartreeFockFieldProj_del(vbb_field)
		CALL SymHartreeFockFieldProj_del(vc_field)
		CALL SymHartreeFockFieldProj_del(vls_field)
		CALL SymHartreeFockFieldProj_del(Delta)

		RETURN
	END SUBROUTINE SymGammaTensorVAP_get_DeltaPhi
	
	!-------------------------------------------------------------------------------!
	!										!
	!   	CALCULATING THE PROJECTED HFB MATRIX: PAIRING DELTA			!
	!										!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGammaTensorVAP_get_DeltaP(HF_Delta_P, CoeffsXY, C_Matrix, NGauge, Delta01)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_Delta_P
		TYPE (ProjectionCoeffs), INTENT(IN) :: CoeffsXY
		INTEGER, INTENT(IN) :: NGauge

		TYPE (SymHartreeFockFieldProj), DIMENSION(:), POINTER :: Delta01
		TYPE (SymHartreeFockFieldProj) :: HF, field1, C_field
		
		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix
				
		COMPLEX, DIMENSION(1:Nmax, 1:Nmax) :: M, MT
		INTEGER :: IndexGauge_tau, ta

		! Calculate the HFB matrix projected field
		
		CALL SymHartreeFockFieldProj_new(HF_Delta_P)
		
		CALL SymHartreeFockFieldProj_new(HF)
		CALL SymHartreeFockFieldProj_new(field1)
		CALL SymHartreeFockFieldProj_new(C_field)
						
		DO IndexGauge_tau = 1, NGauge
		
			!CALL SymGenDensityHFProj_assign_VAP(C_field, C_Matrix(IndexGauge_tau), 0)
			
			CALL SymGenDensityHFProj_assign2(C_field, C_Matrix(IndexGauge_tau))

			DO ta = 0, 1
		
				CALL SymHartreeFockFieldProj_full_product_iso(HF, Delta01(IndexGauge_tau), C_field, ta, ta)
				
				CALL SymHartreeFockFieldProj_scalar_iso(field1, -1.0*CoeffsXY%y_l(IndexGauge_tau, ta), HF, ta)
								
				CALL SymHartreeFockFieldProjIso_add(HF_Delta_P, HF_Delta_P, field1, ta)
				
			END DO
			
		END DO
		
		CALL SymHartreeFockFieldProj_del(HF)
		CALL SymHartreeFockFieldProj_del(field1)
		CALL SymHartreeFockFieldProj_del(C_field)
		
		RETURN
	END SUBROUTINE SymGammaTensorVAP_get_DeltaP

	
	!-------------------------------------------------------------------------------!
	!										!
	!   	CALCULATING THE PROJECTED HFB MATRIX: MEAN-FIELD GAMMA			!
	!										!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGammaTensorVAP_get_GammaP(HF_Gamma_P, CoeffsXY, genden_gauge, C_Matrix, Deriv_Y, NGauge, Gamma, nucleus)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_Gamma_P
		TYPE (ProjectionCoeffs), INTENT(IN) :: CoeffsXY
		TYPE (NucleusType), INTENT(IN) :: nucleus
		INTEGER, INTENT(IN) :: NGauge

		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix, Deriv_Y
		
		TYPE (SymHartreeFockFieldProj), DIMENSION(:, :), POINTER :: GammaP
		TYPE (SymHartreeFockFieldProj), DIMENSION(:, :), POINTER :: Gamma
		
		TYPE (SymHartreeFockFieldProj) :: HF, HF_sum, field1, field2, dY_field, rho_gauge, rho_degen
		
		TYPE (SymDensityProj), DIMENSION(:), POINTER :: density_gauge
		TYPE (SymDensityProj) :: dY
				
		COMPLEX :: Trace_1, Trace_2, Sum
		DOUBLE PRECISION :: AngleGauge, pi
		INTEGER :: IndexGauge_tau, IndexGauge_tau_prime, IndexGauge, N, Z, ta, isoCoupling

		N = Nucleus_get_N(nucleus)
		Z = Nucleus_get_Z(nucleus)

		! Create the HFB matrix gauge-dependent fields
		
		ALLOCATE(GammaP(0:1, 1:NGauge))
		
		pi = 4.0*ATAN(1.0)
		
		DO isoCoupling = 0, 1
		
			DO IndexGauge = 1, NGauge
		
				AngleGauge = pi*REAL(IndexGauge)/REAL(NGauge)
		
				CALL SymHartreeFockFieldProj_new(GammaP(isoCoupling, IndexGauge))
		
				CALL SymGammaTensorVAP_get_GammaProj(genden_gauge, C_Matrix, isoCoupling, AngleGauge, IndexGauge, Gamma, GammaP)
		
			END DO
			
		END DO
		
		! Transform the density into the proper form before to do anything 
		! (required because when "converting" densities, both isospins are treated
		! simultaneously)
		
		ALLOCATE(density_gauge(1:NGauge))
		
		DO IndexGauge_tau = 1, NGauge
			CALL SymDensity_new_GenDensityProj_only_rho(density_gauge(IndexGauge_tau), genden_gauge(IndexGauge_tau))			
		END DO
				
		! Calculate the HFB matrix projected field
		
		CALL SymHartreeFockFieldProj_new(HF_Gamma_P)
		
		CALL SymHartreeFockFieldProj_new(HF)
		CALL SymHartreeFockFieldProj_new(field1)
		CALL SymHartreeFockFieldProj_new(field2)
		CALL SymHartreeFockFieldProj_new(dY_field)
		CALL SymHartreeFockFieldProj_new(rho_gauge)
		CALL SymHartreeFockFieldProj_new(rho_degen)
		
		DO IndexGauge_tau = 1, NGauge
		
			CALL SymGenDensityHFProj_assign_VAP(dY_field, Deriv_Y(IndexGauge_tau), 0)
			CALL SymGenDensityHFProj_assign_VAP(rho_gauge, genden_gauge(IndexGauge_tau)%rho, 0)
			CALL SymGenDensityHFProj_assign_VAP(rho_degen, genden_gauge(IndexGauge_tau)%rho, 1)
			
			DO ta = 0, 1
							
				! First Term
				
				isoCoupling = 0 ! terms that do not couple isospins
				
				Trace_1 = SymHartreeFockFieldProj_product_iso(Gamma(isoCoupling, IndexGauge_tau), rho_degen, ta)
			
				CALL SymHartreeFockFieldProj_scalar_iso(field1, 0.5*Trace_1, dY_field, ta)
				
				CALL SymHartreeFockFieldProj_scalar_iso(field2, CoeffsXY%y_l(IndexGauge_tau, ta), &
										GammaP(isoCoupling, IndexGauge_tau), ta)
				
				CALL SymHartreeFockFieldProjIso_add(HF, field1, field2, ta)
			
				! Loop over isospin 1-ta (tau' in the notes)
			
				isoCoupling = 1 ! terms that DO couple isospins
				
				CALL SymHartreeFockFieldProj_new(HF_sum)
				
				Sum = CMPLX(0,0)
				
				DO IndexGauge_tau_prime = 1, NGauge
				
					Trace_2 = SymHartreeFockFieldProj_product_iso(Gamma(isoCoupling, IndexGauge_tau_prime), rho_degen, ta)

					Sum = Sum + CoeffsXY%y_l(IndexGauge_tau_prime, 1-ta) * Trace_2
					
					CALL SymHartreeFockFieldProj_scalar_iso(field1, CoeffsXY%y_l(IndexGauge_tau_prime, 1-ta), &
										GammaP(isoCoupling, IndexGauge_tau_prime), ta)

					CALL SymHartreeFockFieldProjIso_add(HF_sum, HF_sum, field1, ta)
				
				
				END DO
				
				CALL SymHartreeFockFieldProj_scalar_iso(field1, 0.5*Sum, dY_field, ta)
				
				CALL SymHartreeFockFieldProj_scalar_iso(field2, CoeffsXY%y_l(IndexGauge_tau, ta), HF_sum, ta)
				
				CALL SymHartreeFockFieldProjIso_add(HF_sum, field1, field2, ta)
				
				
				CALL SymHartreeFockFieldProjIso_add(field1, HF, HF_sum, ta)
				CALL SymHartreeFockFieldProjIso_add(HF_Gamma_P, HF_Gamma_P, field1, ta)
				
				CALL SymHartreeFockFieldProj_del(HF_sum)
		
			END DO
			
		END DO
		
		CALL SymHartreeFockFieldProj_del(rho_degen)
		CALL SymHartreeFockFieldProj_del(rho_gauge)
		CALL SymHartreeFockFieldProj_del(HF)
		CALL SymHartreeFockFieldProj_del(field1)
		CALL SymHartreeFockFieldProj_del(field2)
		CALL SymHartreeFockFieldProj_del(dY_field)
		
		DEALLOCATE(density_gauge)
		
		RETURN
	END SUBROUTINE SymGammaTensorVAP_get_GammaP

	!-------------------------------------------------------------------------------!
	!										!
	!   	CALCULATING THE PROJECTED HFB MATRIX: MEAN-FIELD LAMBDA			!
	!										!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGammaTensorVAP_get_LambdaP(HF_Lambda_P, CoeffsXY, genden_gauge, C_Matrix, Deriv_Y, NGauge, Delta01, Delta10, nucleus)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_Lambda_P
		TYPE (ProjectionCoeffs), INTENT(IN) :: CoeffsXY
		TYPE (NucleusType), INTENT(IN) :: nucleus
		INTEGER, INTENT(IN) :: NGauge

		TYPE (SymHartreeFockFieldProj), DIMENSION(:), POINTER :: Delta01, Delta10
		TYPE (SymHartreeFockFieldProj), DIMENSION(:), POINTER :: Lambda
		
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix, Deriv_Y
		
		TYPE (SymHartreeFockFieldProj) :: HF, field1, field2, rho_degen, dY_field
		
		TYPE (SymDensityProj), DIMENSION(:), POINTER :: density_gauge
		TYPE (SymDensityProj) :: dY
				
		COMPLEX :: Trace_1, Trace_2
		DOUBLE PRECISION :: AngleGauge, pi
		INTEGER :: IndexGauge_tau, N, Z, ta

		N = Nucleus_get_N(nucleus)
		Z = Nucleus_get_Z(nucleus)

		! Create the HFB matrix gauge-dependent fields
		
		ALLOCATE(Lambda(1:NGauge))
		
		pi = 4.0*ATAN(1.0)
		
		DO IndexGauge_tau = 1, NGauge
		
			AngleGauge = pi*REAL(IndexGauge_tau)/REAL(NGauge)
		
			CALL SymHartreeFockFieldProj_new(Lambda(IndexGauge_tau))
		
			CALL SymGammaTensorVAP_get_LambdaProj(genden_gauge, C_Matrix, AngleGauge, IndexGauge_tau, Delta01, Delta10, Lambda)
		
		END DO
		
		! Transform the density into the proper form before to do anything 
		! (required because when "converting" densities, both isospins are treated
		! simultaneously)
		
		ALLOCATE(density_gauge(1:NGauge))
		
		DO IndexGauge_tau = 1, NGauge
			CALL SymDensity_new_GenDensityProj01(density_gauge(IndexGauge_tau), genden_gauge(IndexGauge_tau))			
		END DO
		 
		! Calculate the HFB matrix projected field
		
		CALL SymHartreeFockFieldProj_new(HF_Lambda_P)
		
		CALL SymHartreeFockFieldProj_new(field1)
		CALL SymHartreeFockFieldProj_new(field2)
		CALL SymHartreeFockFieldProj_new(dY_field)
		CALL SymHartreeFockFieldProj_new(rho_degen)
		
		DO IndexGauge_tau = 1, NGauge
		
			CALL SymGenDensityHFProj_assign_VAP(dY_field, Deriv_Y(IndexGauge_tau), 0)
			CALL SymGenDensityHFProj_assign_VAP(rho_degen, genden_gauge(IndexGauge_tau)%kap01, 1)
			
			DO ta = 0, 1

				! First Term
				
				CALL SymHartreeFockFieldProj_new(HF)
				
				!Trace_1 = SymHartreeFockFieldProj_product_iso(Delta10(IndexGauge_tau), &
				!					   density_gauge(IndexGauge_tau)%field%kap, ta)
				Trace_1 = SymHartreeFockFieldProj_product_iso(Delta10(IndexGauge_tau), rho_degen, ta)
			
				CALL SymHartreeFockFieldProj_scalar_iso(field1, -0.5*Trace_1, dY_field, ta)
				
				CALL SymHartreeFockFieldProj_scalar_iso(field2, CoeffsXY%y_l(IndexGauge_tau, ta), Lambda(IndexGauge_tau), ta)
				
				CALL SymHartreeFockFieldProjIso_add(HF, field1, field2, ta)

				CALL SymHartreeFockFieldProjIso_add(HF_Lambda_P, HF_Lambda_P, field2, ta)
		
				CALL SymHartreeFockFieldProj_del(HF)
				
			END DO
		END DO
		
		CALL SymHartreeFockFieldProj_del(dY_field)
		CALL SymHartreeFockFieldProj_del(field1)
		CALL SymHartreeFockFieldProj_del(field2)
		CALL SymHartreeFockFieldProj_del(rho_degen)
		
		DEALLOCATE(density_gauge)
		
		RETURN
	END SUBROUTINE SymGammaTensorVAP_get_LambdaP

	!-------------------------------------------------------------------------------!
	!										!
	!  This subroutine calculates:							!
	!										!
	!	   |      d				     |				!
	!     Trace| 0.5*---- Gamma_{tau, tau'}(phi)rho(phi) |				!
	!	   |     d rho				     |				!
	!										!
	!  for both cases tau = tau' and tau different from tau'			!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGammaTensorVAP_get_GammaProj(genden_gauge, C_Matrix, isoCoupling, AngleGauge, IndexGauge, Gamma, GammaP)
		DOUBLE PRECISION, INTENT(IN) :: AngleGauge
		INTEGER, INTENT(IN) :: isoCoupling, IndexGauge

		TYPE (SymHartreeFockFieldProj), DIMENSION(:, :), POINTER :: GammaP
		TYPE (SymHartreeFockFieldProj), DIMENSION(:, :), POINTER :: Gamma
		
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix
		
		TYPE (SymHartreeFockFieldProj) :: HF_tmp_1, HF_tmp_2, HF_out, C_field, rho_field
		TYPE (SymDensityProj) :: density_rho
		
		COMPLEX :: Phase
		INTEGER :: tb, ta
		
		Phase = - 1.0 + CEXP(-CMPLX(0,1)*2.0*AngleGauge)
		
		CALL SymDensity_new_GenDensityProj_only_rho(density_rho, genden_gauge(IndexGauge))

		CALL SymHartreeFockFieldProj_new(HF_tmp_1)
		CALL SymHartreeFockFieldProj_new(HF_tmp_2)
		CALL SymHartreeFockFieldProj_new(HF_out)
		CALL SymHartreeFockFieldProj_new(C_field)
		CALL SymHartreeFockFieldProj_new(rho_field)
		
		CALL SymGenDensityHFProj_assign2(C_field, C_Matrix(IndexGauge))
		!CALL SymGenDensityHFProj_assign_VAP(C_field, C_Matrix(IndexGauge), 1)
		CALL SymGenDensityHFProj_assign_VAP(rho_field, genden_gauge(IndexGauge)%rho, 0)
		
		DO ta = 0, 1
		
						tb = 1-ta
			IF (isoCoupling .EQ. 0) tb = ta
		
			CALL SymHartreeFockFieldProj_full_product_iso(HF_tmp_1, Gamma(isoCoupling, IndexGauge), C_field, ta, tb)
			CALL SymHartreeFockFieldProj_full_product_iso(HF_tmp_2, rho_field, HF_tmp_1, ta, tb)
		
			CALL SymHartreeFockFieldProj_scalar_iso(HF_out, Phase, HF_tmp_2, ta)
		
			CALL SymHartreeFockFieldProjIso_add(GammaP(isoCoupling, IndexGauge), HF_tmp_1, HF_out, ta)
		
		END DO
		
		CALL SymHartreeFockFieldProj_del(rho_field)
		CALL SymHartreeFockFieldProj_del(HF_tmp_1)
		CALL SymHartreeFockFieldProj_del(HF_tmp_2)
		CALL SymHartreeFockFieldProj_del(HF_out)
		CALL SymHartreeFockFieldProj_del(C_field)
		
		RETURN
	END SUBROUTINE SymGammaTensorVAP_get_GammaProj

	!-------------------------------------------------------------------------------!
	!										!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGammaTensorVAP_get_GammaDDProj(genden_gauge, C_Matrix, ta, AngleGauge, IndexGauge, Gamma, GammaP)
		DOUBLE PRECISION, INTENT(IN) :: AngleGauge
		INTEGER, INTENT(IN) :: ta, IndexGauge

		TYPE (SymHartreeFockFieldProj),INTENT(OUT) :: GammaP
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: Gamma
		
		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		
		TYPE (SymHartreeFockFieldProj) :: rho_gauge, HF_tmp_1, HF_tmp_2, HF_out, C_field
		
		COMPLEX :: Phase
		
		Phase = - 1.0 + CEXP(-CMPLX(0,1)*2.0*AngleGauge)
		
		CALL SymHartreeFockFieldProj_new(rho_gauge)
		CALL SymHartreeFockFieldProj_new(C_field)
		
		CALL SymHartreeFockFieldProj_new(HF_tmp_1)
		CALL SymHartreeFockFieldProj_new(HF_tmp_2)
		CALL SymHartreeFockFieldProj_new(HF_out)
		
		
		CALL SymGenDensityHFProj_assign_VAP(C_field, C_Matrix(IndexGauge), 0)
		CALL SymGenDensityHFProj_assign_VAP(rho_gauge, genden_gauge(IndexGauge)%rho, 0)

		CALL SymHartreeFockFieldProj_full_product_iso(HF_tmp_1, Gamma, C_field, ta, ta)
		CALL SymHartreeFockFieldProj_full_product_iso(HF_tmp_2, rho_gauge, HF_tmp_1, ta, ta)
		
		CALL SymHartreeFockFieldProj_scalar_iso(HF_out, Phase, HF_tmp_2, ta)
		
		CALL SymHartreeFockFieldProjIso_add(GammaP, HF_tmp_1, HF_out, ta)
		
		
		CALL SymHartreeFockFieldProj_del(HF_tmp_1)
		CALL SymHartreeFockFieldProj_del(HF_tmp_2)
		CALL SymHartreeFockFieldProj_del(HF_out)
				
		CALL SymHartreeFockFieldProj_del(C_field)
		CALL SymHartreeFockFieldProj_del(rho_gauge)
			
		
		RETURN
	END SUBROUTINE SymGammaTensorVAP_get_GammaDDProj

	
	!-------------------------------------------------------------------------------!
	!										!
	!  This subroutine calculates:							!
	!										!
	!-------------------------------------------------------------------------------!

	SUBROUTINE SymGammaTensorVAP_get_LambdaProj(genden_gauge, C_Matrix, AngleGauge, IndexGauge, Delta01, Delta10, Lambda)
		DOUBLE PRECISION, INTENT(IN) :: AngleGauge
		INTEGER, INTENT(IN) :: IndexGauge

		TYPE (SymHartreeFockFieldProj), DIMENSION(:), POINTER :: Lambda
		TYPE (SymHartreeFockFieldProj), DIMENSION(:), POINTER :: Delta01, Delta10
		
		TYPE (SymGenDensityHFProj), DIMENSION(:), POINTER :: C_Matrix
		TYPE (SymGenDensityGaugeProj), DIMENSION(:), POINTER :: genden_gauge
		
		TYPE (SymHartreeFockFieldProj) :: HF_tmp_1, HF_tmp_2, HF_out, C_field, kap10_field, kap01_field
		TYPE (SymDensityProj) :: density_kap01, density_kap10
		
		COMPLEX :: Phase
		INTEGER :: ta
		
		Phase = 0.5*(1.0 - CEXP(-CMPLX(0,1)*2.0*AngleGauge))
		
		!CALL SymDensity_new_GenDensityProj01(density_kap01, genden_gauge(IndexGauge))
		!CALL SymDensity_new_GenDensityProj10(density_kap10, genden_gauge(IndexGauge))
		
		
		CALL SymHartreeFockFieldProj_new(HF_tmp_1)
		CALL SymHartreeFockFieldProj_new(HF_tmp_2)
		CALL SymHartreeFockFieldProj_new(HF_out)
		CALL SymHartreeFockFieldProj_new(C_field)
		
		CALL SymHartreeFockFieldProj_new(kap10_field)
		CALL SymHartreeFockFieldProj_new(kap01_field)
		
		CALL SymGenDensityHFProj_assign_VAP(C_field, C_Matrix(IndexGauge), 0) 
		CALL SymGenDensityHFProj_assign_VAP(kap10_field, genden_gauge(IndexGauge)%kap10, 0)
		CALL SymGenDensityHFProj_assign_VAP(kap01_field, genden_gauge(IndexGauge)%kap01, 0)

		DO ta = 0, 1
		
			CALL SymHartreeFockFieldProj_full_product_iso(HF_tmp_1, Delta01(IndexGauge), C_field, ta, ta)
			!CALL SymHartreeFockFieldProj_full_product_iso(HF_tmp_2, density_kap10%field%kap, HF_tmp_1, ta, ta)
			CALL SymHartreeFockFieldProj_full_product_iso(HF_tmp_2, kap10_field, HF_tmp_1, ta, ta)
		
			CALL SymHartreeFockFieldProj_full_product_iso(HF_tmp_1, Delta10(IndexGauge), C_field, ta, ta)
			!CALL SymHartreeFockFieldProj_full_product_iso(HF_out, density_kap01%field%kap, HF_tmp_1, ta, ta)
			CALL SymHartreeFockFieldProj_full_product_iso(HF_out, kap01_field, HF_tmp_1, ta, ta)
		
			CALL SymHartreeFockFieldProjIso_add(HF_tmp_1, HF_tmp_2, HF_out, ta)
		
			CALL SymHartreeFockFieldProj_scalar_iso(Lambda(IndexGauge), Phase, HF_tmp_1, ta)

		END DO
		
		CALL SymHartreeFockFieldProj_del(kap10_field)
		CALL SymHartreeFockFieldProj_del(kap01_field)
		
		CALL SymHartreeFockFieldProj_del(HF_tmp_1)
		CALL SymHartreeFockFieldProj_del(HF_tmp_2)
		CALL SymHartreeFockFieldProj_del(HF_out)
		CALL SymHartreeFockFieldProj_del(C_field)
		
		RETURN
	END SUBROUTINE SymGammaTensorVAP_get_LambdaProj

	
	

END MODULE selfc_proj
