MODULE wave

	USE nucleus
	USE symfield
	USE genden

	IMPLICIT NONE

	TYPE WaveFunction
		TYPE (NucleusType) nucleus
		TYPE (GenDensity) genden
		INTEGER num
!		INTEGER dimensionUV
	END TYPE

CONTAINS

	SUBROUTINE WaveFunction_new(wave, N, Z)
		TYPE (WaveFunction), INTENT(INOUT) :: wave
		INTEGER, INTENT(IN) :: N, Z

		CALL Nucleus_new(wave%nucleus, N, Z)
		CALL GenDensity_new(wave%genden)
		wave%num = 0
		RETURN
	END SUBROUTINE WaveFunction_new

	SUBROUTINE WaveFunction_new_Nucleus(wave, nuc)
		TYPE (WaveFunction), INTENT(INOUT) :: wave
		TYPE (NucleusType), INTENT(IN) :: nuc

		CALL Nucleus_new_Nucleus(wave%nucleus, nuc)
		RETURN
	END SUBROUTINE WaveFunction_new_Nucleus

	SUBROUTINE WaveFunction_write(wave, n)
		TYPE (WaveFunction), INTENT(INOUT) :: wave
		INTEGER, INTENT(IN) :: n

		INTEGER file_error
		INTEGER, PARAMETER :: file_id = 6

		PRINT *
		PRINT *, "Escribiendo resultados en:", wave%nucleus%filename
! filename = filename + ".0"
!		OPEN(file_id, FILE=wave%nucleus%filename, ACTION="WRITE", IOSTAT=file_error)
!		IF (file_error .NE. 0) STOP "ERROR: Imposible escribir fichero"

!		CALL Nucleus_write(wave%nucleus, file_id)
!		CALL SymHartreeFockField_write(wave%genden%rho, file_id)
!		CALL SymHartreeFockField_write(wave%genden%kap, file_id)
!		CLOSE(file_id)
		RETURN
	END SUBROUTINE WaveFunction_write

END MODULE wave
