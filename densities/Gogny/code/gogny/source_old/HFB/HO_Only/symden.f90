 MODULE symden

	USE input
	USE math
	USE global
	USE symd3t
	USE nucleus
	USE symfield

	IMPLICIT NONE

	!---------------------------------------------------------------!
	!    SymDensity is a type for objects associated to:		!
	!	- a nucleus 						!
	!       - a couple of mean-field and pairing field Gamma and 	!
	!         Delta							!
	!---------------------------------------------------------------!

	TYPE SymDensity
		TYPE (NucleusType) nucleus
		TYPE (SymHartreeFockBogolField) field
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE SymDensity_assign
	END INTERFACE

 CONTAINS

	SUBROUTINE SymDensity_new(density, N, Z)
		! Recibe como parametros de entrada el registro de densidad
		! y el numero de protones y electrones
		TYPE (SymDensity), INTENT(INOUT) :: density
		INTEGER, INTENT(IN) :: N, Z
		DOUBLE PRECISION :: b

		CALL SymHartreeFockBogolField_new(density%field)
		CALL Nucleus_new(density%nucleus, N, Z, b)
		
		CALL SymDensity_read(density)
		
		RETURN
	END SUBROUTINE SymDensity_new

	SUBROUTINE SymDensity_new_SymDensity(density_out, density_in)
		TYPE (SymDensity), INTENT(INOUT) :: density_out
		TYPE (SymDensity), INTENT(IN) :: density_in

		CALL SymHartreeFockBogolField_new(density_out%field)
		CALL Nucleus_new_Nucleus(density_out%nucleus, density_in%nucleus)
		
		density_out%field%rho = density_in%field%rho
		density_out%field%kap = density_in%field%kap
		
		RETURN
	END SUBROUTINE SymDensity_new_SymDensity

	SUBROUTINE SymDensity_initialize(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		INTEGER :: ta, la, na, i, lower_bound
		INTEGER :: nosc, fin_apa, ini_apa, num, den, lup
		DOUBLE PRECISION :: vv, vu

		! Calculo de la densidad inicial
		DO ta = 0, 1
		
			! Initialization of density matrices
			density%field%rho%p(ta) = DBLE(0.0)
			density%field%rho%a(ta) = DBLE(0.0)
			density%field%kap%p(ta) = DBLE(0.0)
			density%field%kap%a(ta) = DBLE(0.0)

			! Actual number of particles for isospin ta
			num = density%nucleus%np(ta)
			i = 1
			IF (N_0 .LE. 1) THEN
				IF (num .GE. MagicNumber(N_0)) STOP "Abortado: SymDensity_initialize (num >= MagicNumber(N_0))"
			END IF

			DO WHILE (MagicNumber(i) .LE. num)
				i = i + 1
			END DO
			fin_apa = i

			IF (ta .EQ. 1) THEN
				PRINT "(/A,I2,A,I3)", "Hay apareamiento de los ", density%nucleus%np(ta), " neutrones hasta el estado mas bajo de la capa ", i
			ELSE
				PRINT "(/A,I2,A,I3)", "Hay apareamiento de los ", density%nucleus%np(ta), " protones hasta el estado mas bajo de la capa ", i
			END IF

 			den = MagicNumber(i)
			vv = DBLE(num) / DBLE(den)
   			vu = SQRT(vv) * SQRT(1.0 - vv)
			DO nosc = i - 1, 0, -1
				! Convenio u = (-) ^ l * |u|
				vu = vu * PAR(nosc)
                                lup = MIN(nosc, N_0)
				DO la = lup, 0, -2
					na = ((nosc - la) / 2) + 1
					density%field%rho%p(ta)%d3tensor(la)%d2(na, na) = 2 * ((2 * la) + 1) * vv
					density%field%kap%p(ta)%d3tensor(la)%d2(na, na) = 2 * ((2 * la) + 1) * vu
				END DO
			END DO

			la = MIN(fin_apa, N_0)
			na = 1
			density%field%rho%p(ta)%d3tensor(la)%d2(na, na) = 2 * (la + 1) * vv
			density%field%rho%a(ta)%d3tensor(la)%d2(na, na) = vv
			density%field%kap%p(ta)%d3tensor(la)%d2(na, na) = 2 * (la + 1) * vu
			density%field%kap%a(ta)%d3tensor(la)%d2(na, na) = vu
					
			density%nucleus%actual_np(ta) = SymD3Tensor_trace(density%field%rho%p(ta))
		END DO
		
		RETURN
	END SUBROUTINE SymDensity_initialize

	SUBROUTINE SymDensity_read(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		IF(SymHartreeFockBogolField_read(density%field, density%nucleus%filename)) THEN
			density%nucleus%actual_np(0) = SymD3Tensor_trace(density%field%rho%p(0))
			density%nucleus%actual_np(1) = SymD3Tensor_trace(density%field%rho%p(1))
		ELSE
			CALL SymDensity_initialize(density)
		END IF
		RETURN
	END SUBROUTINE SymDensity_read

	SUBROUTINE SymDensity_save(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		CALL SymHartreeFockBogolField_write(density%field, density%nucleus%filename)
		
		RETURN
	END SUBROUTINE SymDensity_save
	
	SUBROUTINE SymDensity_assign(density_out, density_in)
		TYPE (SymDensity), INTENT(INOUT) :: density_out
		TYPE (SymDensity), INTENT(IN) :: density_in

		density_out%field = density_in%field
		RETURN
	END SUBROUTINE SymDensity_assign

	SUBROUTINE SymDensity_store_actual_R2(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		DOUBLE PRECISION b2
		INTEGER ta

		b2 = Nucleus_get_b(density%nucleus) ** 2
		
		DO ta=0, 1
			density%nucleus%actual_R2(ta) = R2Field * density%field%rho%p(ta) / density%nucleus%np(ta)
		END DO
		
		density%nucleus%actual_R2(2) = &
			 ((density%nucleus%np(0) * density%nucleus%actual_R2(0))  &
			+ (density%nucleus%np(1) * density%nucleus%actual_R2(1))) &
			/ (density%nucleus%np(0) + density%nucleus%np(1))
			
		RETURN
	END SUBROUTINE SymDensity_store_actual_R2

	SUBROUTINE SymDensity_shuffle(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		INTEGER ta, la, nosc
		INTEGER i, fin_apa, ini_apa, num, den, na
		DOUBLE PRECISION vv, vu

		DO ta = 0, 1
		
			density%field%kap%p(ta) = DBLE(0.0)
			density%field%kap%a(ta) = DBLE(0.0)
			num = density%nucleus%np(ta)
			
			IF (N_0 .LE. 8) THEN
				IF (num .GE. MagicNumber(N_0)) STOP "Numero de particulas no valido"
			END IF
			i = 1
			DO WHILE (MagicNumber(i) .LE. num)
				i = i + 1
			END DO
			fin_apa = i

			den = MagicNumber(i)
			vv = DBLE(num) / den
			vu = SQRT(vv) * SQRT(1.0 - vv)
			
			DO nosc = i - 1, 0, -1
				vu = vu * PAR(nosc) ! convenio u = (-)^l * |u|
				DO la = nosc, 0, -2
					na = ((nosc - la) / 2) + 1
					density%field%kap%p(ta)%d3tensor(la)%d2(na, na) = 2.0 * ((2.0 * la) + 1.0) * vu
				END DO
			END DO
			la = fin_apa
			na = 1

			density%field%kap%p(ta)%d3tensor(la)%d2(na, na) = 2.0 * (la + 1.0) * vu
			density%field%kap%a(ta)%d3tensor(la)%d2(na, na) = vu
		END DO
		RETURN
	END SUBROUTINE SymDensity_shuffle

	SUBROUTINE SymDensity_show_SpatialDistribution(density, ta)
		TYPE (SymDensity), INTENT(INOUT) :: density
		INTEGER, INTENT(IN) :: ta

!TODO
		RETURN
	END SUBROUTINE SymDensity_show_SpatialDistribution

	SUBROUTINE SymDensity_show_ParticleDensity(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

!TODO Donde esta?
		RETURN
	END SUBROUTINE SymDensity_show_ParticleDensity

	SUBROUTINE SymDensity_del(density)
		TYPE (SymDensity), INTENT(INOUT) :: density

		CALL Nucleus_del(density%nucleus)
		CALL SymHartreeFockBogolField_del(density%field)
		
		RETURN
	END SUBROUTINE SymDensity_del

END MODULE symden
