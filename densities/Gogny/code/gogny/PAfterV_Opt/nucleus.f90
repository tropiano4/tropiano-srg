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

	CHARACTER(5), DIMENSION(0:102) :: Nucleus_specie
	DATA Nucleus_specie / "n", &
		"H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", &
		"Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",        &
		"Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",       &
		"Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",       &
		"Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",       &
		"Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", &
		"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", &
		"Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", &
		"Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", &
		"Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu",             &
		"Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No"/
 CONTAINS

	!---------------------------------------------------------------------------------------!
	! This subroutine initializes an object of type nucleus. It predefines quantities such	!
	! as Z and N, the Fermi levels, the oscillator length, etc. It also defines blocked 	!
	! level as function of the number of particles.						!
	!---------------------------------------------------------------------------------------!

	SUBROUTINE Nucleus_new(nuc, N, Z, b)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, INTENT(IN) :: N, Z
		DOUBLE PRECISION, INTENT(IN), OPTIONAL:: b

		CHARACTER(LEN = 32) :: specie, N_0_str, A_str
		INTEGER, DIMENSION(1:200) :: ShellModel_L, ShellModel_J

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
		nuc%lambda_np = (/-35.0, -0.5     /)
		nuc%actual_R2 = (/ 0.0,  0.0, 0.0/)
		nuc%lambda_R2 = (/ 0.0,  0.0     /)

		IF (PRESENT(b)) THEN
			nuc%func%x = b
		ELSE
			nuc%func%x = Nucleus_get_InitialOscillatorLength(nuc)
		END IF

		! Get the shell model spectrum to define "blocked" levels
		CALL Nucleus_ShellModel(ShellModel_L, ShellModel_J)
		
		! Defining blocked levels for neutrons or protons (depending)
		CALL Nucleus_set_NeutronBlocking(nuc, ShellModel_L, ShellModel_J)
		CALL Nucleus_set_ProtonBlocking(nuc, ShellModel_L, ShellModel_J)

		! Name of the file that contains the density 
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

	!---------------------------------------------------------------------------------------!
	! This subroutine copies the content of one "nucleus" into another			!
	!---------------------------------------------------------------------------------------!

	SUBROUTINE Nucleus_new_Nucleus(n1, n2)
		TYPE (NucleusType), INTENT(INOUT) :: n1
		TYPE (NucleusType), INTENT(IN) :: n2

		CALL R1R1Function_new(n1%func)
		CALL Nucleus_copy(n1, n2)
		
		RETURN
	END SUBROUTINE Nucleus_new_Nucleus

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

	!-----------------------------------------------------------------------!
	!  In the case of a harmonic oscillator basis, we predefine the oscil-	!
	!  lator length	as: b = 1.01*A**(1/6)					!
	!-----------------------------------------------------------------------!

	FUNCTION Nucleus_get_InitialOscillatorLength(nuc)
		DOUBLE PRECISION Nucleus_get_InitialOscillatorLength
		TYPE (NucleusType), INTENT(IN) :: nuc

		Nucleus_get_InitialOscillatorLength = 1.01 * (DBLE(nuc%np(0) + nuc%np(1)) ** (1.0 / 6.0))
		
		RETURN
	END FUNCTION Nucleus_get_InitialOscillatorLength

	!-----------------------------------------------------------------------!
	!  We define which level is blocked for the protons and fill in all the	!
	!  relevant quantities of that level					!
	!-----------------------------------------------------------------------!

	SUBROUTINE Nucleus_set_ProtonBlocking(nuc, ShellModel_L, ShellModel_J)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, DIMENSION(1:200), INTENT(IN) :: ShellModel_L, ShellModel_J

		INTEGER :: Z, N, A
		INTEGER :: lp, jp, ip, mu0p

		IF (nuc%num .EQ. 0) THEN
			nuc%is_blocking(0) = .FALSE.
			RETURN
		END IF

		Z = nuc%np(0)
		N = nuc%np(1)
		A = nuc%np(0) + nuc%np(1)

		nuc%is_blocking(0) = IS_ODD(Z)
		
		IF (nuc%is_blocking(0)) THEN
		
			lp = ShellModel_L(Z)
			jp = ShellModel_J(Z)

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

	!-----------------------------------------------------------------------!
	!  We define which level is blocked for the neutrons and fill in all the!
	!  relevant quantities of that level					!
	!-----------------------------------------------------------------------!

	SUBROUTINE Nucleus_set_NeutronBlocking(nuc, ShellModel_L, ShellModel_J)
		TYPE (NucleusType), INTENT(INOUT) :: nuc
		INTEGER, DIMENSION(1:200), INTENT(IN) :: ShellModel_L, ShellModel_J

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
		
			ln = ShellModel_L(N)
			jn = ShellModel_J(N)

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

	!-----------------------------------------------------------------------------------------------!
	! This subroutine defines the shell-model orbitals that we will use in the blocking		!
	! For each orbital, we give the orbital angular momentum L and the total angular momentum J	!
	!-----------------------------------------------------------------------------------------------!

	SUBROUTINE Nucleus_ShellModel(ShellModel_L, ShellModel_J)
		INTEGER, DIMENSION(1:200), INTENT(OUT) :: ShellModel_L, ShellModel_J

		INTEGER ln, jn, i 

		! 1s1/2
		DO i=1,2
			ShellModel_L(i) = 0
			ShellModel_J(i) = 1
		END DO

		! 1p3/2
		DO i=3,6
			ShellModel_L(i) = 1
			ShellModel_J(i) = 3
		END DO

		! 1p1/2
		DO i=7,8
			ShellModel_L(i) = 1
			ShellModel_J(i) = 1
		END DO

		! 1d5/2
		DO i=9,14
			ShellModel_L(i) = 2
			ShellModel_J(i) = 5
		END DO

		! 2s1/2
		DO i=15,16
			ShellModel_L(i) = 0
			ShellModel_J(i) = 1
		END DO

		! 1d3/2
		DO i=17,20
			ShellModel_L(i) = 2
			ShellModel_J(i) = 3
		END DO

		! 1f7/2
		DO i=21,28
			ShellModel_L(i) = 3
			ShellModel_J(i) = 7
		END DO

		! 2p3/2
		DO i=29,32
			ShellModel_L(i) = 1
			ShellModel_J(i) = 3
		END DO

		! 1f5/2
		DO i=33,38
			ShellModel_L(i) = 3
			ShellModel_J(i) = 5
		END DO

		! 2p1/2
		DO i=39,40
			ShellModel_L(i) = 1
			ShellModel_J(i) = 1
		END DO

		! 1g9/2
		DO i=41,50
			ShellModel_L(i) = 4
			ShellModel_J(i) = 9
		END DO

		! 1g7/2
		DO i=51,58
			ShellModel_L(i) = 4
			ShellModel_J(i) = 7
		END DO

		! 2d5/2
		DO i=59,64
			ShellModel_L(i) = 2
			ShellModel_J(i) = 5
		END DO

		! 2d3/2
		DO i=65,68
			ShellModel_L(i) = 2
			ShellModel_J(i) = 3
		END DO

		! 3s1/2
		DO i=69,70
			ShellModel_L(i) = 0
			ShellModel_J(i) = 1
		END DO

		! 1h11/2
		DO i=71,82
			ShellModel_L(i) = 5
			ShellModel_J(i) = 11
		END DO

		! 1h9/2
		DO i=83,92
			ShellModel_L(i) = 5
			ShellModel_J(i) = 9
		END DO

		! 2f7/2
		DO i=93,100
			ShellModel_L(i) = 3
			ShellModel_J(i) = 7
		END DO

		RETURN
	END SUBROUTINE Nucleus_ShellModel

	!-------------------------------------------------------------------------------!
	! Various utilities follow. Taking as input an object of type "Nucleus", it:	!
	!   - returns the mass number (Nucleus_get_A)					!
	!   - returns the oscillator length (Nucleus_get_b)				!
	!   - returns the neutron number (Nucleus_get_N)				!
	!   - returns the proton number (Nucleus_get_Z)					!
	!   - returns the current rms mass radius (Nucleus_get_actual_R2)		!
	!   - returns the square of the current neutron r.m.s. radius (Nucleus_get_R2n)	!
	!   - returns the square of the current proton r.m.s. radius (Nucleus_get_R2p)	!
	!   - sets the neutron number (Nucleus_set_N)					!
	!   - sets the proton number (Nucleus_set_Z)					!
	!   - sets the oscillator length (Nucleus_set_b)				!
	!-------------------------------------------------------------------------------!

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

	! Devuelve el numero de neutrones del nucleo
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

	!-------------------------------------------------------------------------------!
	!  Subroutine that displays all what is contained in a given "Nucleus" object	!
	!-------------------------------------------------------------------------------!
	
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
