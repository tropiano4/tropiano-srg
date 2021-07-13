! Input parameters module
 MODULE input

	IMPLICIT NONE

	INTEGER :: Npoint

	DOUBLE PRECISION, ALLOCATABLE :: Radius(:)
	DOUBLE PRECISION, ALLOCATABLE :: DensityNeutron(:),DensityProton(:)
	DOUBLE PRECISION, ALLOCATABLE :: WaveFunUProt(:),WaveFunVProt(:),WaveFunUNeut(:),WaveFunVNeut(:)

 CONTAINS

        SUBROUTINE ReadWaveFunction()

		INTEGER, PARAMETER :: file_unit = 13
		INTEGER :: file_error, stat, Nline, i

		DOUBLE PRECISION :: dummy1,dummy2,dummy3,dummy4,dummy5

		OPEN(file_unit, FILE="WavesQP.dat", ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) STOP "Impossible to open: WavesQP.dat"

		Nline = 0
		DO
			Nline = Nline + 1
			READ(file_unit, *, IOSTAT=stat) dummy1,dummy2,dummy3,dummy4,dummy5
			IF (stat.NE.0) EXIT
		END DO

		Npoint = Nline-1  ! Number of points minus the first one (=0)

		REWIND (file_unit)

		IF (.NOT. ALLOCATED(Radius)) 		ALLOCATE(Radius(1:Npoint))
		IF (.NOT. ALLOCATED(WaveFunUNeut)) 	ALLOCATE(WaveFunUNeut(1:Npoint))
		IF (.NOT. ALLOCATED(WaveFunVNeut)) 	ALLOCATE(WaveFunVNeut(1:Npoint))
		IF (.NOT. ALLOCATED(WaveFunUProt)) 	ALLOCATE(WaveFunUProt(1:Npoint))
		IF (.NOT. ALLOCATED(WaveFunVProt)) 	ALLOCATE(WaveFunVProt(1:Npoint))

		READ(file_unit, *, IOSTAT=stat) (Radius(i),WaveFunVProt(i),WaveFunUProt(i),WaveFunVNeut(i),WaveFunUNeut(i),i=1,Npoint)

		!IF (ABS(WaveFunUProt(i)).LE.1.E-9) DensityProton(i)=0.0

		RETURN
	END SUBROUTINE ReadWaveFunction

        SUBROUTINE ReadDensity()

		INTEGER, PARAMETER :: file_unit = 13
		INTEGER :: file_error, stat, Nline, i

		DOUBLE PRECISION :: dummy1,dummy2,dummy3

		OPEN(file_unit, FILE="density.txt", ACTION="READ", IOSTAT=file_error)
		IF (file_error .NE. 0) STOP "Impossible to open: density.txt"

		Nline = 0
		DO
			Nline = Nline + 1
			READ(file_unit, *, IOSTAT=stat) dummy1,dummy2,dummy3
			IF (stat.NE.0) EXIT
		END DO

		Npoint = Nline - 3 ! Number of points minus the first one (=0)

		REWIND (file_unit)

		IF (.NOT. ALLOCATED(Radius)) 		ALLOCATE(Radius(1:Npoint))
		IF (.NOT. ALLOCATED(DensityProton)) 	ALLOCATE(DensityProton(1:Npoint))
		IF (.NOT. ALLOCATED(DensityNeutron)) 	ALLOCATE(DensityNeutron(1:Npoint))

		READ(file_unit, *, IOSTAT=stat) dummy1,dummy2,dummy3
		READ(file_unit, *, IOSTAT=stat) (Radius(i),DensityProton(i),DensityNeutron(i),i=1,Npoint)

		IF (ABS(DensityProton(i)).LE.1.E-9) DensityProton(i)=0.0

		RETURN
	END SUBROUTINE ReadDensity

	SUBROUTINE PrintOut(Xvector, Yvector)
		DOUBLE PRECISION, INTENT(IN) :: Xvector(:), Yvector(:)

		INTEGER :: i

		DO i = 1, Npoint
			WRITE(6,'(2f20.14)') Xvector(i),Yvector(i)
		END DO

		RETURN

	END SUBROUTINE PrintOut

END MODULE input
