PROGRAM main

        !---------------------------------------------------------------------!
        !    								      !
        !                      INCLUDING MODULES                              !
        !    								      !
        !---------------------------------------------------------------------!
	
	USE input    ! Parametros de entrada
	USE global
	USE nucleus  ! Nucleo
	USE symden
	USE diagmeth ! Metodo de diagonalizacion

        !---------------------------------------------------------------------!
        !    								      !
        !                 DECLARATION OF VARIABLES                            !
        !    								      !
        !---------------------------------------------------------------------!
	
	IMPLICIT NONE
	
        TYPE (SymDensity) density
	TYPE (DiagonalizationMethod) diagonal

	INTEGER cycles_in, cycles_out, cycles_rate
	DOUBLE PRECISION min_b
		
        !---------------------------------------------------------------------!
        !    								      !
        !                 BEGINNING OF THE PROGRAM                            !
        !    								      !
        !---------------------------------------------------------------------!

	CALL SYSTEM_CLOCK(cycles_in)

	CALL Global_new

	CALL SymDensity_new(density, neutrons, protons)
	CALL Nucleus_set_b(density%nucleus, b_0)

	CALL DiagonalizationMethod_new(diagonal, density)
	CALL DiagonalizationMethod_goto_SelfConsistency(diagonal, Convergence)

	CALL Global_del

	CALL SYSTEM_CLOCK(cycles_out, cycles_rate)
	PRINT "(/A)", "Tiempo consumido:"
	PRINT "(A,F12.5)", "Segundos: ", DBLE(cycles_out - cycles_in) / cycles_rate
	PRINT "(A,I5)", "Ciclos de reloj: ", cycles_out - cycles_in
	
	STOP
	
END PROGRAM main
