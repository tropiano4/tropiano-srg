	! This is only used to define a type "OneDimSolve" that will 
	! be used in the determination of the particle number. This file 
	! is not used as module, but directly in an INCLUDE statement.

	TYPE OneDimSolve
		TYPE (R1R1Function), POINTER :: func
		DOUBLE PRECISION AccX ! Accuracy in X direction
		DOUBLE PRECISION AccY ! Accuracy in Y direction
		INTEGER lim ! Maximum number of iterations
		DOUBLE PRECISION, DIMENSION(0:2) :: x, y ! Trial values of X and Y
		INTEGER L, C, U, Last ! Locations of trial values
		INTEGER vpol, hpol ! Polarities
		DOUBLE PRECISION YY ! Target value
		LOGICAL Finish ! .TRUE. if LookAt finds conv.
		LOGICAL Captured ! .TRUE. when target surrounded
		TYPE (DiagonalizationMethod), POINTER :: diagonal
	END TYPE
