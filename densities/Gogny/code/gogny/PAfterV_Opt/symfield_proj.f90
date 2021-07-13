 MODULE symfield_proj

	USE input
	USE global
	USE symd3t_proj
	USE symtalm

	IMPLICIT NONE

	!-------------------------------------------------------------------------------!
	!  SymHartreeFockFieldProj :Made of 4 tensors (2 for protons, 2 for neutrons)	!
	!    - p refers to the local part of the HF field: Gamma_local			!
	!    - a refers to exchange part of the HF field: Gamma_exchange		!
	!-------------------------------------------------------------------------------!

	TYPE SymHartreeFockFieldProj
		TYPE (SymD3Tensor), DIMENSION(0:1) :: p
		TYPE (SymD3Tensor), DIMENSION(0:1) :: a
		DOUBLE PRECISION, DIMENSION(0:1) :: GaugeAngle
	END TYPE

	!-------------------------------------------------------------------------------!
	!  SymHartreeFockBogoFieldProj: Made of 2 SymHartreeFockFieldProj tensors 	!
	!     - rho refers to the p-h mean-field density: Gamma (local and exchange)	!
	!     - kap refers to p-p pairing field: Delta	 (local and exchange)		!
	!-------------------------------------------------------------------------------!

	TYPE SymHartreeFockBogolFieldProj
		TYPE (SymHartreeFockFieldProj) rho
		TYPE (SymHartreeFockFieldProj) kap
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE &
			SymHartreeFockFieldProj_assign, &
			SymHartreeFockBogolFieldProj_assign
	END INTERFACE

	INTERFACE OPERATOR(*)
		MODULE PROCEDURE SymHartreeFockFieldProj_product2
	END INTERFACE

 CONTAINS

	SUBROUTINE SymHartreeFockBogolFieldProj_new(HFB)
		TYPE (SymHartreeFockBogolFieldProj), INTENT(INOUT) :: HFB

		CALL SymHartreeFockFieldProj_new(HFB%rho)
		CALL SymHartreeFockFieldProj_new(HFB%kap)

		RETURN
	END SUBROUTINE SymHartreeFockBogolFieldProj_new

	SUBROUTINE SymHartreeFockBogolFieldProj_setGauge(HFB, Gauge, ta)
		TYPE (SymHartreeFockBogolFieldProj), INTENT(INOUT) :: HFB
		DOUBLE PRECISION, INTENT(IN) :: Gauge
		INTEGER, INTENT(IN) ::  ta

		CALL SymHartreeFockFieldProj_setGauge(HFB%rho, Gauge, ta)
		CALL SymHartreeFockFieldProj_setGauge(HFB%kap, Gauge, ta)

		RETURN
	END SUBROUTINE SymHartreeFockBogolFieldProj_setGauge

	FUNCTION SymHartreeFockBogolFieldProj_read(HFB, filename)
		LOGICAL SymHartreeFockBogolFieldProj_read

		TYPE (SymHartreeFockBogolFieldProj), INTENT(INOUT) :: HFB

		CHARACTER(*), INTENT(IN) :: filename

		INTEGER, PARAMETER :: file_desc = 16 ! Descriptor de fichero
		INTEGER :: file_error, new_N_0, new_Lmax, new_Nmax, new_Basis, new_CompHO

		OPEN (file_desc, FILE=TRIM(filename), STATUS="OLD", ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			SymHartreeFockBogolFieldProj_read = .FALSE.
			RETURN
		END IF

		READ (file_desc, *) new_N_0, new_Lmax, new_Nmax, new_Basis, new_CompHO

		IF (new_Basis .NE. Basis) THEN
			CLOSE (file_desc)
			SymHartreeFockBogolFieldProj_read = .FALSE.
			RETURN
		ELSE
			IF (Basis .EQ. 1) THEN ! Case HO analytical
				IF (new_N_0 .NE. N_0) THEN
					CLOSE (file_desc)
					SymHartreeFockBogolFieldProj_read = .FALSE.
					RETURN
				END IF
			ELSE ! Case WS or HO numerical
				IF (new_compHO .NE. CompHO) THEN
					CLOSE (file_desc)
					SymHartreeFockBogolFieldProj_read = .FALSE.
					RETURN
				ELSE
					IF (compHO .EQ. 1 .AND. new_N_0 .NE. N_0) THEN
						CLOSE (file_desc)
						SymHartreeFockBogolFieldProj_read = .FALSE.
						RETURN
					END IF
					IF (compHO .EQ. 0 .AND. (new_Lmax .NE. Lmax .OR. new_Nmax .NE. Nmax)) THEN
						CLOSE (file_desc)
						SymHartreeFockBogolFieldProj_read = .FALSE.
						RETURN
					END IF
				END IF
			END IF
		END IF

		IF ((SymHartreeFockFieldProj_read(HFB%rho, file_desc, file_error)) .AND. &
		    (SymHartreeFockFieldProj_read(HFB%kap, file_desc, file_error))) THEN
			SymHartreeFockBogolFieldProj_read = .TRUE.
		ELSE
			SymHartreeFockBogolFieldProj_read = .FALSE.
		END IF
		CLOSE (file_desc)

		RETURN
	END FUNCTION SymHartreeFockBogolFieldProj_read

	! Almacena los datos de los campos rho y kap en un fichero
	SUBROUTINE SymHartreeFockBogolFieldProj_write(HBF, filename)
		TYPE (SymHartreeFockBogolFieldProj), INTENT(IN) :: HBF
		CHARACTER(*), INTENT(IN) :: filename

		INTEGER, PARAMETER :: file_desc = 16 ! Descriptor de fichero
		INTEGER file_error

		OPEN (file_desc, FILE=TRIM(filename), ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) RETURN

		! El primer valor almacenado en el fichero es el numero de capas: N_0
		WRITE (file_desc, FMT="(5I3)", IOSTAT=file_error) N_0, Lmax, Nmax, Basis, CompHO

		! A continuacion se almacenan rho y kap
		CALL SymHartreeFockFieldProj_write(HBF%rho, file_desc, file_error)
		CALL SymHartreeFockFieldProj_write(HBF%kap, file_desc, file_error)

		CLOSE (file_desc)
		RETURN
	END SUBROUTINE SymHartreeFockBogolFieldProj_write

	SUBROUTINE SymHartreeFockBogolFieldProj_assign(HBF1, HBF2)
		TYPE (SymHartreeFockBogolFieldProj), INTENT(INOUT) :: HBF1
		TYPE (SymHartreeFockBogolFieldProj), INTENT(IN) :: HBF2

		HBF1%rho = HBF2%rho
		HBF1%kap = HBF2%kap

		RETURN
	END SUBROUTINE SymHartreeFockBogolFieldProj_assign

	SUBROUTINE SymHartreeFockBogolFieldProj_del(HFB)
		TYPE (SymHartreeFockBogolFieldProj), INTENT(INOUT) :: HFB

		CALL SymHartreeFockFieldProj_del(HFB%rho)
		CALL SymHartreeFockFieldProj_del(HFB%kap)

		RETURN
	END SUBROUTINE SymHartreeFockBogolFieldProj_del

	SUBROUTINE SymHartreeFockFieldProj_new(HF)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF

		INTEGER ta

		DO ta = 0, 1

			CALL SymD3Tensor_new(HF%p(ta))
			CALL SymD3Tensor_new(HF%a(ta))

			HF%GaugeAngle(ta) = DBLE(-999.9)

		END DO

		RETURN
	END SUBROUTINE SymHartreeFockFieldProj_new

	FUNCTION SymHartreeFockFieldProj_read(HF, file_desc, file_error)
		LOGICAL SymHartreeFockFieldProj_read

		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF

		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		INTEGER ta, length

		DO ta = 0, 1
			IF (.NOT. SymD3Tensor_read(HF%p(ta), file_desc, file_error)) THEN
				SymHartreeFockFieldProj_read = .FALSE.
				RETURN
			END IF
			IF (.NOT. SymD3Tensor_read(HF%a(ta), file_desc, file_error)) THEN
				SymHartreeFockFieldProj_read = .FALSE.
				RETURN
			END IF
		END DO

		SymHartreeFockFieldProj_read = .TRUE.

		RETURN
	END FUNCTION SymHartreeFockFieldProj_read

	SUBROUTINE SymHartreeFockFieldProj_write(HF, file_desc, file_error)
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF

		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		INTEGER ta

		DO ta = 0, 1
			CALL SymD3Tensor_write(HF%p(ta), file_desc, file_error)
			CALL SymD3Tensor_write(HF%a(ta), file_desc, file_error)
		END DO

		RETURN
	END SUBROUTINE SymHartreeFockFieldProj_write

	SUBROUTINE SymHartreeFockFieldProj_assign(HF_out, HF_in)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF_in

		INTEGER ta

		DO ta = 0, 1

			HF_out%p(ta) = HF_in%p(ta)
			HF_out%a(ta) = HF_in%a(ta)

			HF_out%GaugeAngle(ta) = HF_in%GaugeAngle(ta)

		END DO

		RETURN
	END SUBROUTINE SymHartreeFockFieldProj_assign

	SUBROUTINE SymHartreeFockFieldProj_setGauge(HF_out, Gauge, ta)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		DOUBLE PRECISION, INTENT(IN) :: Gauge
		INTEGER, INTENT(IN) ::  ta

		HF_out%GaugeAngle(ta) = Gauge

		RETURN
	END SUBROUTINE SymHartreeFockFieldProj_setGauge

	SUBROUTINE SymHartreeFockFieldProj_add(HF_out, HF1_in, HF2_in)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF1_in, HF2_in

		INTEGER :: ta

		CALL SymD3Tensor_add(HF_out%p(0), HF1_in%p(0), HF2_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), HF1_in%p(1), HF2_in%p(1))

		CALL SymD3Tensor_add(HF_out%a(0), HF1_in%a(0), HF2_in%a(0))
		CALL SymD3Tensor_add(HF_out%a(1), HF1_in%a(1), HF2_in%a(1))

		DO ta = 0, 1
			!IF (HF1_in%GaugeAngle(ta) .NE. HF2_in%GaugeAngle(ta)) THEN
			!	WRITE(*,'("Gauge angle 1 = ",F12.4)') HF1_in%GaugeAngle(ta)
			!	WRITE(*,'("Gauge angle 2 = ",F12.4)') HF2_in%GaugeAngle(ta)
			!	STOP "Error in SymHartreeFockFieldProj_add of module symfieldProj.f90 - Different Gauge Angles"
			!ELSE
				HF_out%GaugeAngle(ta) = HF1_in%GaugeAngle(ta)
			!END IF
		END DO

		RETURN
	END SUBROUTINE SymHartreeFockFieldProj_add

	SUBROUTINE SymHartreeFockFieldProj_add_SymD3Tensor(HF_out, t_in, HF_in)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out

		TYPE (SymD3Tensor), INTENT(IN) :: t_in
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF_in

		CALL SymD3Tensor_add(HF_out%p(0), t_in, HF_in%p(0))
		CALL SymD3Tensor_add(HF_out%p(1), t_in, HF_in%p(1))

		HF_out%a(0) = HF_in%a(0)
		HF_out%a(1) = HF_in%a(1)

		HF_out%GaugeAngle(0) = HF_in%GaugeAngle(0)
		HF_out%GaugeAngle(1) = HF_in%GaugeAngle(1)

		RETURN
	END SUBROUTINE SymHartreeFockFieldProj_add_SymD3Tensor

	SUBROUTINE SymHartreeFockFieldProj_product(HT_out, R1, HF_in)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HT_out

		DOUBLE PRECISION, INTENT(IN) :: R1
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF_in

		CALL SymD3Tensor_product(HT_out%p(0), R1, HF_in%p(0))
		CALL SymD3Tensor_product(HT_out%p(1), R1, HF_in%p(1))

		CALL SymD3Tensor_product(HT_out%a(0), R1, HF_in%a(0))
		CALL SymD3Tensor_product(HT_out%a(1), R1, HF_in%a(1))

		HT_out%GaugeAngle(0) = HF_in%GaugeAngle(0)
		HT_out%GaugeAngle(1) = HF_in%GaugeAngle(1)

		RETURN
	END SUBROUTINE SymHartreeFockFieldProj_product

	! Used to make the contraction of a field of given isospin with
	! a density of the same isospin, then summing over both isospins.
	! The result is a real number which does NOT depend on the
	! isospin.

	FUNCTION SymHartreeFockFieldProj_product2(HF1_in, HF2_in)
		DOUBLE PRECISION SymHartreeFockFieldProj_product2
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF1_in, HF2_in

		DOUBLE PRECISION, DIMENSION(0:1) :: sum2
		INTEGER ta

		DO ta = 0, 1
			sum2(ta) = (HF1_in%p(ta) * HF2_in%p(ta)) &
			         + (HF1_in%a(ta) * HF2_in%a(ta))
		END DO

		SymHartreeFockFieldProj_product2 = sum2(0) + sum2(1)

		RETURN
	END FUNCTION SymHartreeFockFieldProj_product2

	! Used to make the contraction of a field of given isospin with
	! a density of the same isospin. The result is a real number
	! (which depends on the isospin).

	FUNCTION SymHartreeFockFieldProj_product_iso(HF1_in, HF2_in, ta)
		DOUBLE PRECISION SymHartreeFockFieldProj_product_iso
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF1_in, HF2_in

		DOUBLE PRECISION :: sum2
		INTEGER ta

		sum2 = (HF1_in%p(ta) * HF2_in%p(ta)) &
		     + (HF1_in%a(ta) * HF2_in%a(ta))

		SymHartreeFockFieldProj_product_iso = sum2

		RETURN
	END FUNCTION SymHartreeFockFieldProj_product_iso

	FUNCTION SymHartreeFockFieldProj_distance(HF1_in, HF2_in)
		DOUBLE PRECISION SymHartreeFockFieldProj_distance
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF1_in, HF2_in

		SymHartreeFockFieldProj_distance = MAX( &
			SymD3Tensor_distance(HF1_in%p(0), HF2_in%p(0)), &
			SymD3Tensor_distance(HF1_in%p(1), HF2_in%p(1)))

		RETURN
	END FUNCTION SymHartreeFockFieldProj_distance

	FUNCTION SymHartreeFockFieldProj_ChargeDensity(HF, r)
		DOUBLE PRECISION SymHartreeFockFieldProj_ChargeDensity
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF
		DOUBLE PRECISION, INTENT(IN) :: r

		DOUBLE PRECISION x, d1, xk, sum1, sumlb
		INTEGER k, lb, nb, nd, nbmax, p2max, p2

		x = r * r
		xk = DBLE(1.0)
		sum1 = DBLE(0.0)

		DO k = 0, Lmax

			sumlb = DBLE(0.0)

			DO lb = 0, k
									nbmax = MIN(Nmax, NmaxOfL(lb))
				IF (basis .EQ. 1 .OR. CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

				DO nb = 1, nbmax
					DO nd = 1, nbmax
						p2max = nb + nd - 2
						p2 = k - lb
						IF (p2 .GT. p2max) CYCLE
						d1 = SymCoefficientB_get(nb - 1, lb, nd - 1, lb, p2)
						sumlb = sumlb + (d1 * HF%p(PROTON)%d3tensor(lb)%d2(nd, nb))
					END DO
				END DO
			END DO
			sum1 = sum1 + (xk * sumlb)
			xk = xk * x
		END DO
		SymHartreeFockFieldProj_ChargeDensity = EXP(-x) * sum1
		RETURN
	END FUNCTION SymHartreeFockFieldProj_ChargeDensity

	SUBROUTINE SymHartreeFockFieldProj_del(HF)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF

		INTEGER ta

		DO ta = 0, 1
			CALL SymD3Tensor_del(HF%p(ta))
			CALL SymD3Tensor_del(HF%a(ta))
		END DO

		RETURN
	END SUBROUTINE SymHartreeFockFieldProj_del

END MODULE symfield_proj
