PROGRAM main

	USE input    ! Parametros de entrada
	USE global
	USE nucleus  ! Nucleo
	USE symden   ! Matriz simetrica de densidad
	USE diagmeth ! Metodo de diagonalizacion
	USE wave

	IMPLICIT NONE

	TYPE (SymDensity) density
	TYPE (DiagonalizationMethod) diagonal

	INTEGER cycles_in, cycles_out, cycles_rate
	DOUBLE PRECISION min_b

	CALL SYSTEM_CLOCK(cycles_in)

	CALL Global_new

	CALL SymDensity_new(density, protons, neutrons) ! N y Z
	CALL Nucleus_set_b(density%nucleus, b_0)
	CALL SymDensity_read(density)

	CALL DiagonalizationMethod_new(diagonal, density)
	CALL DiagonalizationMethod_goto_SelfConsistency(diagonal, DBLE(1.0e-4))

	CALL Global_del

	CALL SYSTEM_CLOCK(cycles_out, cycles_rate)
	PRINT "(/A)", "Tiempo consumido:"
	PRINT "(A,EN)", "Segundos: ", DBLE(cycles_out - cycles_in) / cycles_rate
	PRINT "(A,I)", "Ciclos de reloj: ", cycles_out - cycles_in
	STOP
END PROGRAM main
