PROGRAM main

	USE input    ! Parametros de entrada
	USE global
	USE nucleus  ! Nucleo
	USE symden
	USE diagmeth ! Metodo de diagonalizacion

	IMPLICIT NONE

        TYPE (SymDensity) :: density
	TYPE (DiagonalizationMethod) :: diagonal

	INTEGER :: cycles_in, cycles_out, cycles_rate

	! Initialize clock
	CALL SYSTEM_CLOCK(cycles_in)

	! Define important quantities and computes one-body kinetic energy and
	! Gauss-Laguerre quadratures (global.f90)
	CALL Global_new

	! Define a new set of densities and a new nucleus
	CALL SymDensity_new(density, neutrons, protons)
	CALL Nucleus_set_b(density%nucleus, b_0)

	! Self-consistent loop
	CALL DiagonalizationMethod_new(diagonal, density)
	CALL DiagonalizationMethod_goto_SelfConsistency(diagonal, Convergence)

	! Nothing happens here
	CALL Global_del

	! Stops the clock and print times
	CALL SYSTEM_CLOCK(cycles_out, cycles_rate)
	WRITE(6,'("Elapsed Time")')
	WRITE(6,'("Seconds .............:",F20.6)') DBLE(cycles_out - cycles_in) / cycles_rate
	WRITE(6,'("Clock cycles ........:",I12)') cycles_out - cycles_in

	STOP

END PROGRAM main
