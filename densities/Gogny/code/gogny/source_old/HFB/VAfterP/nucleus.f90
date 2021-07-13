MODULE nucleus

	USE math
	USE r1r1

	IMPLICIT NONE

	CHARACTER(LEN = 10), DIMENSION(0:1) :: Nucleus_type
	DATA Nucleus_type / "protones", "neutrones" /

	TYPE NucleusType
		TYPE (R1R1Function) func
		INTEGER, DIMENSION(0:2) :: np
		DOUBLE PRECISION, DIMENSION(0:1) :: actual_np, lambda_np, lambda_R2
		DOUBLE PRECISION, DIMENSION(0:2) :: actual_R2
		DOUBLE PRECISION eHFB, R2, norma
		DOUBLE PRECISION pairing ! Apareamiento
		LOGICAL is_blocking(0:1)
		INTEGER num
		INTEGER, DIMENSION(0:1) :: ia, la, ja, mu0
		CHARACTER(LEN = 64) :: filename
	END TYPE

	INTERFACE ASSIGNMENT(=)
		MODULE PROCEDURE Nucleus_copy
	END INTERFACE

	INCLUDE "name_new.f90"
	!INCLUDE "binding.f90"
	!INCLUDE "atomic.f90"
	!INCLUDE "mass.f90"

 CONTAINS

	SUBROUTINE Nucleus_new(nuc, N, Z, b)
		! Recibe como parámetros de entrada
		! el registro de densidad y el número de protones y electrones
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, INTENT(IN) :: N, Z
		DOUBLE PRECISION, INTENT(IN), OPTIONAL:: b

		CHARACTER(LEN = 32) :: specie, N_0_str, A_str

		CALL R1R1Function_new(nuc%func)

		nuc%eHFB    = 0.0
		nuc%pairing = 0.0
		nuc%R2      = 0.0
		nuc%norma   = 1.0
		nuc%num     = 1
		nuc%np(0)   = Z
		nuc%np(1)   = N
		nuc%np(2)   = N + Z
		nuc%actual_np = (/ 0.0,  0.0     /)
		nuc%lambda_np = (/-10.0, -10.0     /)
		nuc%actual_R2 = (/ 0.0,  0.0, 0.0/)
		nuc%lambda_R2 = (/ 0.0,  0.0     /)

		IF (PRESENT(b)) THEN
			nuc%func%x = b
		ELSE
			nuc%func%x = Nucleus_get_InitialOscillatorLength(nuc)
		END IF

		CALL Nucleus_set_NeutronBlocking(nuc)
		CALL Nucleus_set_ProtonBlocking(nuc)

		! Nombre de fichero
		IF (nuc%np(0) .LE. 102) THEN
			specie = Nucleus_specie(nuc%np(0))
		ELSE
			CALL Int2Char(specie, nuc%np(0))
		END IF
		CALL Int2Char(A_str, nuc%np(2))
		CALL Int2Char(N_0_str, N_0)
		nuc%filename = "data/" // TRIM(specie) // TRIM(A_str) // "_" // TRIM(N_0_str)
		RETURN
	END SUBROUTINE Nucleus_new

	SUBROUTINE Nucleus_new_Nucleus(n1, n2)
		TYPE (NucleusType), INTENT(INOUT) :: n1
		TYPE (NucleusType), INTENT(IN) :: n2

		CALL R1R1Function_new(n1%func)
		CALL Nucleus_copy(n1, n2)
		RETURN
	END SUBROUTINE Nucleus_new_Nucleus

	FUNCTION Nucleus_read(nucleus_out, file_desc, file_error)
		LOGICAL Nucleus_read
		TYPE (NucleusType), INTENT(INOUT) :: nucleus_out
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		DOUBLE PRECISION x

		Nucleus_read = .FALSE.
		READ (file_desc, FMT="(E)", IOSTAT=file_error) x
		IF (file_error .NE. 0) RETURN
		CALL R1R1Function_set(nucleus_out%func, x)

		READ (file_desc, FMT="(E)", IOSTAT=file_error) nucleus_out%eHFB
		IF (file_error .NE. 0) RETURN

		READ (file_desc, FMT="(E)", IOSTAT=file_error) nucleus_out%R2
		IF (file_error .NE. 0) RETURN

		READ (file_desc, FMT="(E)", IOSTAT=file_error) nucleus_out%pairing
		IF (file_error .NE. 0) RETURN

		READ (file_desc, FMT="(E,E)", IOSTAT=file_error) nucleus_out%actual_R2(1), nucleus_out%actual_R2(0)
		IF (file_error .NE. 0) RETURN

		READ (file_desc, FMT="(E,E)", IOSTAT=file_error) nucleus_out%lambda_R2(1), nucleus_out%lambda_R2(0)
		IF (file_error .NE. 0) RETURN

		READ (file_desc, FMT="(E,E)", IOSTAT=file_error) nucleus_out%actual_np(1), nucleus_out%actual_np(0)
		IF (file_error .NE. 0) RETURN

		READ (file_desc, FMT="(E,E)", IOSTAT=file_error) nucleus_out%lambda_np(1), nucleus_out%lambda_np(0)
		IF (file_error .NE. 0) RETURN

		READ (file_desc, FMT="(E)", IOSTAT=file_error) nucleus_out%norma
		IF (file_error .NE. 0) RETURN

		Nucleus_read = .TRUE.
		RETURN
	END FUNCTION Nucleus_read

	SUBROUTINE Nucleus_write(nucleus_in, file_desc, file_error)
		TYPE (NucleusType), INTENT(IN) :: nucleus_in
		INTEGER, INTENT(IN) :: file_desc
		INTEGER, INTENT(INOUT) :: file_error

		WRITE (file_desc, FMT="(E)", IOSTAT=file_error) nucleus_in%func%x
		WRITE (file_desc, FMT="(E)", IOSTAT=file_error) nucleus_in%eHFB
		WRITE (file_desc, FMT="(E)", IOSTAT=file_error) nucleus_in%R2
		WRITE (file_desc, FMT="(E)", IOSTAT=file_error) nucleus_in%pairing
		WRITE (file_desc, FMT="(E,E)", IOSTAT=file_error) nucleus_in%actual_R2(1), nucleus_in%actual_R2(0)
		WRITE (file_desc, FMT="(E,E)", IOSTAT=file_error) nucleus_in%lambda_R2(1), nucleus_in%lambda_R2(0)
		WRITE (file_desc, FMT="(E,E)", IOSTAT=file_error) nucleus_in%actual_np(1), nucleus_in%actual_np(0)
		WRITE (file_desc, FMT="(E,E)", IOSTAT=file_error) nucleus_in%lambda_np(1), nucleus_in%lambda_np(0)
		WRITE (file_desc, FMT="(E)", IOSTAT=file_error) nucleus_in%norma
		RETURN
	END SUBROUTINE Nucleus_write

	!-----------------------------------------------------------------------!
	!  Copy the data of nucleus n2 into nucleus n1. Both n1 and n2 have 	!
	!  the derived type NucleusType						!
	!-----------------------------------------------------------------------!

	SUBROUTINE Nucleus_copy(n1, n2)
		TYPE (NucleusType), INTENT(INOUT) :: n1
		TYPE (NucleusType), INTENT(IN) :: n2

		INTEGER ta

		n1%eHFB = n2%eHFB
		n1%pairing = n2%pairing
		n1%R2 = n2%R2
		n1%num = n2%num
		DO ta = 0, 1
			n1%np(ta) = n2%np(ta)
			n1%actual_np(ta) = n2%actual_np(ta)
			n1%lambda_np(ta) = n2%lambda_np(ta)
			n1%actual_R2(ta) = n2%actual_R2(ta)
			n1%lambda_R2(ta) = n2%lambda_R2(ta)
			n1%is_blocking(ta) = n2%is_blocking(ta)
			IF (n1%is_blocking(ta) .EQV. .TRUE.) THEN
				n1%ia(ta) = n2%ia(ta)
				n1%la(ta) = n2%la(ta)
				n1%ja(ta) = n2%ja(ta)
				n1%mu0(ta) = n2%mu0(ta)      
			END IF
		END DO
		n1%np(2) = n1%np(0) + n1%np(1)
		n1%actual_R2(2) = ((n1%np(0) * n1%actual_R2(0)) + (n1%np(1) * n1%actual_R2(1))) / (n1%np(0) + n1%np(1))
  
		CALL Nucleus_set_b(n1, Nucleus_get_b(n2))
		n1%norma = n2%norma
		n1%filename = n2%filename
		RETURN
	END SUBROUTINE Nucleus_copy

	FUNCTION Nucleus_get_InitialOscillatorLength(nuc)
		DOUBLE PRECISION Nucleus_get_InitialOscillatorLength
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_InitialOscillatorLength = 1.01 * (DBLE(nuc%np(0) + nuc%np(1)) ** (1.0 / 6.0))
		RETURN
	END FUNCTION Nucleus_get_InitialOscillatorLength

	SUBROUTINE Nucleus_set_ProtonBlocking(nuc)
		TYPE (NucleusType), INTENT(INOUT) :: nuc

		INTEGER Z, N, A
		INTEGER lp, jp, ip, mu0p

		IF (nuc%num .EQ. 0) THEN
			nuc%is_blocking(0) = .FALSE.
			RETURN
		END IF

		Z = nuc%np(0)
		N = nuc%np(1)
		A = nuc%np(0) + nuc%np(1)

		nuc%is_blocking(0) = IS_ODD(Z)
		
		IF (nuc%is_blocking(0)) THEN
			IF (Z .EQ. 3) THEN
				lp = 1
				jp = 3 ! Li
			ELSE IF ((Z .GT. 8) .AND. (Z .LT. 14)) THEN
				lp = 2
				jp = 5
			ELSE IF (Z .EQ. 15) THEN
				lp = 0
				jp = 1 ! Capa sd
			ELSE IF ((Z .GT. 16) .AND. (Z .LT. 20)) THEN
				lp = 2
				jp = 3
			ELSE IF (Z .EQ. 41) THEN
				lp = 4
				jp = 9 ! Nb
			ELSE IF ((Z .EQ. 23) .OR. (Z .EQ. 27)) THEN
				lp = 3
				jp = 7 ! V y Co
			ELSE
				WRITE(*,'("Odd nucleus but no Proton Blocking defined for N = ",I4," !!!")') N
				STOP 'From Nucleus_set_ProtonBlocking (module nucleus.f90)'				
			END IF

			ip = lp * 2
			IF (jp .EQ. ((2 * lp) - 1)) THEN
				ip = ip - 1
			END IF
			mu0p = 1

			! Se bloquea el estado mas proximo al nivel de Fermi
			WRITE(*,'("Proton Blocking: l = ",I2," j = ",I2,"/2 - mu = ",i2)') lp,jp,mu0p

			nuc%la(0) = lp
			nuc%ja(0) = jp
			nuc%ia(0) = ip
			nuc%mu0(0) = mu0p
		END IF
		RETURN
	END SUBROUTINE Nucleus_set_ProtonBlocking

	SUBROUTINE Nucleus_set_NeutronBlocking(nuc)
		TYPE (NucleusType), INTENT(INOUT) :: nuc

		INTEGER Z, N, A
		INTEGER ln, jn, in, mu0n

		IF (nuc%num .EQ. 0) THEN
			nuc%is_blocking(1) = .FALSE.
			RETURN
		END IF

		Z = nuc%np(0)
		N = nuc%np(1)
		A = nuc%np(0) + nuc%np(1)

		nuc%is_blocking(1) = IS_ODD(N)

		IF (nuc%is_blocking(1)) THEN
			IF ((N .GT. 2) .AND. (N .LT. 6)) THEN
				ln = 1
				jn = 3 ! 1p3/2 state
			ELSE IF (N .EQ. 7) THEN
				ln = 1
				jn = 1 ! 1p1/2 state
			ELSE IF ((N .GT. 8) .AND. (N .LT. 14)) THEN
				ln = 2
				jn = 5 ! 1d5/2 state
			ELSE IF (N .EQ. 15) THEN
				ln = 0
				jn = 1 ! 2s1/2 state
			ELSE IF ((N .GT. 16) .AND. (N .LT. 20)) THEN
				ln = 2
				jn = 3
			ELSE IF ((N .GT. 20) .AND. (N .LT. 28)) THEN
				ln = 3
				jn = 7
			ELSE IF (Z .EQ. 18) THEN
				! Isotopos del Argon
				IF (A .EQ. 51) THEN
					ln = 0
					jn = 1
				ELSE IF ((A .EQ. 53) .OR. (A .EQ. 55)) THEN
					ln = 2
					jn = 3
				ELSE IF ((A .EQ. 57) .OR. (A .EQ. 59)) THEN
					ln = 3
					jn = 7
				ELSE
					WRITE(*,'("Impossible to define Neutron Blocking for Z=18 !!!")')
					STOP 'From Nucleus_set_NeutronBlocking (module nucleus.f90)'
				END IF
			ELSE IF (Z .EQ. 20) THEN
				! Isotopos del Calcio
				IF ((A .EQ. 57) .OR. (A .EQ. 59)) THEN
					ln = 2
					jn = 3
				ELSE IF ((A .GE. 61) .AND. (A .LE. 67)) THEN
					ln = 3
					jn = 7
				ELSE IF ((A .GE. 69) .AND. (A .LE. 73)) THEN
					ln = 1
					jn = 3
				ELSE
					WRITE(*,'("Impossible to define Neutron Blocking for Z=20 !!!")')
					STOP 'From Nucleus_set_NeutronBlocking (module nucleus.f90)'
				END IF
			ELSE IF ((Z .EQ. 38) .AND. (A .EQ. 87)) THEN
				ln = 4
				jn = 9 ! Sr
			ELSE IF ((Z .EQ. 22) .AND. (A .EQ. 49)) THEN
				ln = 3
				jn = 7 ! Ti
			ELSE IF (Z .EQ. 50) THEN ! Sn
				IF ((A .GE. 113) .AND. (A .LE. 119)) THEN
					ln = 0
					jn = 1
				ELSE IF ((A .GE. 121) .AND. (A .LE. 125)) THEN
					ln = 5
					jn = 11
				ELSE
					WRITE(*,'("Impossible to define Neutron Blocking for Z=50 !!!")')
					STOP 'From Nucleus_set_NeutronBlocking (module nucleus.f90)'
				END IF
			ELSE
				WRITE(*,'("Odd nucleus but no Neutron Blocking defined for Z = ",I4," !!!")') Z
				STOP 'From Nucleus_set_NeutronBlocking (module nucleus.f90)'				
			END IF

			in = ln * 2
			IF (jn .EQ. ((2 * ln) - 1)) THEN
				in = in - 1
			END IF
			mu0n = 1

			! Se bloquea el estado mas proximo al nivel de Fermi
			WRITE(*,'("Neutron Blocking: l = ",I2," j = ",I2,"/2 - mu = ",i2)') ln,jn,mu0n

			nuc%la(1) = ln
			nuc%ja(1) = jn
			nuc%ia(1) = in
			nuc%mu0(1) = mu0n
		END IF
		RETURN
	END SUBROUTINE Nucleus_set_NeutronBlocking

	SUBROUTINE Nucleus_show_ExperimentalData(nuc)
		TYPE (NucleusType), INTENT(IN) :: nuc

		CHARACTER(128) str

		IF (nuc%num .EQ. 0) RETURN

		!PRINT "(/A,A)", "Datos experimentales para: ", Nucleus_name(nuc%num)
!		PRINT *, "Binding Energy:", Expm(Nucleus_BindingEnergy(0, nuc%num), Nucleus_BindingEnergy(1, nuc%num), "eV", "k")
!		PRINT *, "Atomic Mass:", Expm(Nucleus_AtomicMass(0, nuc%num), Nucleus_AtomicMass(1, nuc%num), "amu", "micro")
!		PRINT *, "Mass Excess:", Expm(Nucleus_MassExcess(0, nuc%num), Nucleus_MassExcess(1, nuc%num), "eV", "k")
		!PRINT "(A16,E,A,E,A)", "Binding Energy:", Nucleus_BindingEnergy(0, nuc%num), " +/-", Nucleus_BindingEnergy(1, nuc%num), " eV"
		!PRINT "(A16,E,A,E,A)", "Atomic Mass:", Nucleus_AtomicMass(0, nuc%num), " +/-", Nucleus_AtomicMass(1, nuc%num), " amu"
		!PRINT "(A16,E,A,E,A)", "Mass Excess:", Nucleus_MassExcess(0, nuc%num), " +/-", Nucleus_MassExcess(1, nuc%num), " eV"
		RETURN
	END SUBROUTINE Nucleus_show_ExperimentalData

	FUNCTION IS_ODD(n)
		LOGICAL IS_ODD
		INTEGER, INTENT(IN) :: n

		IS_ODD = BTEST(n, 0)
		RETURN
	END FUNCTION IS_ODD

	FUNCTION Nucleus_get_A(nuc)
		INTEGER Nucleus_get_A
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_A = nuc%np(0) + nuc%np(1)
		RETURN
	END FUNCTION Nucleus_get_A

	FUNCTION Nucleus_get_b(nuc)
		DOUBLE PRECISION Nucleus_get_b
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_b = nuc%func%x
		RETURN
	END FUNCTION Nucleus_get_b

	! Devuelve el número de neutrones del núcleo
	FUNCTION Nucleus_get_N(nuc)
		INTEGER Nucleus_get_N
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_N = nuc%np(1)
		RETURN
	END FUNCTION Nucleus_get_N

	! Devuelve el número de protones del núcleo
	FUNCTION Nucleus_get_Z(nuc)
		INTEGER Nucleus_get_Z
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_Z = nuc%np(0)
		RETURN
	END FUNCTION Nucleus_get_Z

	FUNCTION Nucleus_get_actual_R2(nuc)
		DOUBLE PRECISION Nucleus_get_actual_R2
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_actual_R2 = ((nuc%np(0) * nuc%actual_R2(0))  &
		                       + (nuc%np(1) * nuc%actual_R2(1))) &
		                       / DBLE(nuc%np(0) + nuc%np(1))
		RETURN
	END FUNCTION Nucleus_get_actual_R2

	FUNCTION Nucleus_get_R2n(nuc)
		DOUBLE PRECISION Nucleus_get_R2n
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_R2n = nuc%actual_R2(1)
		RETURN
	END FUNCTION Nucleus_get_R2n

	FUNCTION Nucleus_get_R2p(nuc)
		DOUBLE PRECISION Nucleus_get_R2p
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_R2p = nuc%actual_R2(0)
		RETURN
	END FUNCTION Nucleus_get_R2p

	SUBROUTINE Nucleus_set_N(nuc, N)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, INTENT(IN) :: N

		nuc%np(1) = N
		nuc%np(2) = nuc%np(0) + nuc%np(1)
		RETURN
	END SUBROUTINE Nucleus_set_N

	SUBROUTINE Nucleus_set_Z(nuc, Z)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, INTENT(IN) :: Z

		nuc%np(0) = Z
		nuc%np(2) = nuc%np(0) + nuc%np(1)
		RETURN
	END SUBROUTINE Nucleus_set_Z

	SUBROUTINE Nucleus_set_b(nuc, b)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		DOUBLE PRECISION, INTENT(IN) :: b

		CALL R1R1Function_set(nuc%func, b)
		RETURN
	END SUBROUTINE Nucleus_set_b

	! Muestra los valores del núcleo
	SUBROUTINE Nucleus_show_Status(nuc)
		TYPE (NucleusType), INTENT(IN) :: nuc

		INTEGER i

		PRINT *
		PRINT "(A)", "Nucleo:"
		PRINT "(70A1)", ("-", i = 1, 70) ! Muestra una línea de separación

		PRINT "(A30,F15.5)", "Longitud del oscinador (b):", nuc%func%x
		PRINT "(A30,I15,I15)", "Numero de particulas (np):", nuc%np(0), nuc%np(1)
		PRINT "(A30,F15.5,F15.5)", "(actual_np):", nuc%actual_np(0), nuc%actual_np(1)
		PRINT "(A30,F15.5,F15.5)", "(lambda_np):", nuc%lambda_np(0), nuc%lambda_np(1)
		PRINT "(A30,F15.5,F15.5,F15.5)", "(actual_R2):", nuc%actual_R2(0), nuc%actual_R2(1), nuc%actual_R2(2)
		PRINT "(A30,E15.5,E15.5)", "(lambda_R2):", nuc%lambda_R2(0), nuc%lambda_R2(1)
		PRINT "(A30,F15.5)", "Energia total (eHFB):", nuc%eHFB
		PRINT "(A30,E15.5)", "Apareamiento:", nuc%pairing
		PRINT "(A30,E15.5)", "R2:", nuc%R2
		PRINT "(A30,F15.5)", "Norma:", nuc%norma

		PRINT "(70A1)", ("-", i = 1, 70)

		IF (nuc%is_blocking(0)) THEN
			PRINT "(A,I5,A,I5,A,I5)", "Hay bloqueo de neutrones en: la=", nuc%la(0), "ja=", nuc%ja(0), "mu0=", nuc%mu0(0)
		END IF
		IF (nuc%is_blocking(1)) THEN
			PRINT "(A,I5,A,I5,A,I5)", "Hay bloqueo de neutrones en: la=", nuc%la(1), "ja=", nuc%ja(1), "mu0=", nuc%mu0(1)
		END IF
	END SUBROUTINE Nucleus_show_Status

	SUBROUTINE Nucleus_del(nuc)
		TYPE (NucleusType), INTENT(INOUT) :: nuc

		CALL R1R1Function_del(nuc%func)
	END SUBROUTINE Nucleus_del

END MODULE nucleus
