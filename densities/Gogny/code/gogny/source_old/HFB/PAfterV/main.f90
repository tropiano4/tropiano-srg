PROGRAM main

        !---------------------------------------------------------------------!
        !    								      !
        !                      INCLUDING MODULES                              !
        !    								      !
        !---------------------------------------------------------------------!
	
	USE input    ! Parametros de entrada
	USE global
	USE nucleus  ! Nucleo
	USE symden_proj
	USE diagmeth ! Metodo de diagonalizacion

        !---------------------------------------------------------------------!
        !    								      !
        !                 DECLARATION OF VARIABLES                            !
        !    								      !
        !---------------------------------------------------------------------!
	
	IMPLICIT NONE
	
        TYPE (SymDensityProj) :: density
	TYPE (DiagonalizationMethod) :: diagonal

	INTEGER cycles_in, cycles_out, cycles_rate
	DOUBLE PRECISION min_b
		
        !---------------------------------------------------------------------!
        !    								      !
        !                 BEGINNING OF THE PROGRAM                            !
        !    								      !
        !---------------------------------------------------------------------!

	CALL SYSTEM_CLOCK(cycles_in)

	CALL Global_new

	CALL SymDensityProj_new(density, neutrons, protons)
	CALL Nucleus_set_b(density%nucleus, b_0)

	CALL DiagonalizationMethod_new(diagonal, density)
	CALL DiagonalizationMethod_goto_SelfConsistency(diagonal, Convergence)

	CALL Global_del

	CALL SYSTEM_CLOCK(cycles_out, cycles_rate)
	PRINT "(/A)", "Tiempo consumido:"
	PRINT "(A,EN)", "Segundos: ", DBLE(cycles_out - cycles_in) / cycles_rate
	PRINT "(A,I)", "Ciclos de reloj: ", cycles_out - cycles_in
	
	STOP
	
END PROGRAM main
