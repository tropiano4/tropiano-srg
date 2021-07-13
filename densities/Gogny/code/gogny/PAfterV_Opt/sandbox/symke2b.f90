!------------------------------------------------------------------!
!								   !
!  CALCULATION OF THE KINETIC ENERGY TERM OF THE GOGNY FORCE       !
!								   !
!------------------------------------------------------------------!

 MODULE symke2b

	USE input
	USE global
	USE symd3t_proj
	USE symfield_proj

	IMPLICIT NONE

	TYPE SymKineticEnergy2Body
		TYPE (SymD3Tensor_SymD3Tensor) :: v11, v22
	END TYPE

	TYPE SymEk2pp
		TYPE (SymD3Tensor_SymD3Tensor) :: v11_pair, v22_pair
	END TYPE

	CHARACTER(64), PRIVATE :: filename
        INTEGER, PRIVATE :: file_desc = 16

	PRIVATE SymKineticEnergy2Body_nablaHO

 CONTAINS

	SUBROUTINE SymKineticEnergy2Body_new(vEkCMph, vEkCMpp)
		TYPE (SymKineticEnergy2Body), INTENT(INOUT) :: vEkCMph
		TYPE (SymEk2pp), INTENT(INOUT) :: vEkCMpp

		!  Create and initializes the tensors
		CALL SymD3Tensor_SymD3Tensor_new(vEkCMph%v11)
		CALL SymD3Tensor_SymD3Tensor_new(vEkCMph%v22)
		CALL SymD3Tensor_SymD3Tensor_new(vEkCMpp%v11_pair)
		CALL SymD3Tensor_SymD3Tensor_new(vEkCMpp%v22_pair)

                SELECT CASE (Basis)

		CASE (1)

			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/vEkCM", N_0, "ph_HO.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/vEkCM", N_0, "ph_HO.txt"
			END IF

		CASE (2)

			IF (N_0 < 10) THEN
				WRITE(filename, "(A,I1,A)") "data/vEkCM", N_0, "ph_WS.txt"
			ELSE
				WRITE(filename, "(A,I2,A)") "data/vEkCM", N_0, "ph_WS.txt"
			END IF

		END SELECT

	END SUBROUTINE SymKineticEnergy2Body_new

        !---------------------------------------------------------------------------------------!
	!											!
	!  Sum over na, nb, nc, nd and la, lb to calculate the matrix elements of the type:	!
        !											!
	!                   < na la, nb lb | T | nd lb, nc la >					!
        !											!
	!  where T is the kinetic energy. Since T is a one-body operator, this can actually	!
	!  be separated into:									!
	!                        < na la | T | nd lb >< nb lb | T | nc la >			!
	!											!
	!  WARNING: FOR L, THE SUMMATION IS NOT DONE UP TO A GIVEN L_MAX, BUT OVER N_0		!
	!           FOR N, N_MAX IS EXPRESSED AS FUNCTION OF N_0 TOO				!
        !											!
        !---------------------------------------------------------------------------------------!

	SUBROUTINE SymKineticEnergy2Body_calculate(vEkCMph, vEkCMpp)
		TYPE (SymKineticEnergy2Body), INTENT(INOUT) :: vEkCMph
		TYPE (SymEk2pp), INTENT(INOUT) :: vEkCMpp

		INTEGER :: la, lla, na, namax, nc, lb, llb, nb, nbmax, nd
		DOUBLE PRECISION :: cuad, d1

		INTEGER, PARAMETER :: file_desc = 16
		INTEGER :: file_error

		OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)
		IF (file_error .NE. 0) THEN
			PRINT *, "*** Attention: Impossible to write the results in ", filename
		END IF

                SELECT CASE (Basis)

		CASE (1)

			DO la = 0, Lmax

				lla = (2 * la) + 1 ! Degeneracy of a shell with l_a
				namax = ((N_0 - la) / 2) + 1 ! Maximum n value for the "bra"

				DO na = 1, namax
					DO nc = 1, na

			        		DO lb = 0, la

							llb = (2 * lb) + 1 ! Degeneracy of a shell with l_b
							nbmax = ((N_0 - lb) / 2) + 1 ! Maximum n value for the "ket"

							cuad = CUAD2(la, lb, 1) ! Defined in Module math.f90

							DO nb = 1, nbmax
								DO nd = 1, nb

									! Product: < na la || NABLA || nd la >< nb lb || NABLA || nc lb >
									d1 = SymKineticEnergy2Body_nablaHO(na - 1, la, nd - 1, lb) * &
									SymKineticEnergy2Body_nablaHO(nb - 1, lb, nc - 1, la)

									! Anti-symmetrization of the matrix elements by permutation:
									!          < na la || NABLA || nb la >< nd lb || NABLA || nc lb >
									IF (nb .NE. nd) THEN
										d1 = d1 + (SymKineticEnergy2Body_nablaHO(na - 1, la, nb - 1, lb) * &
												SymKineticEnergy2Body_nablaHO(nd - 1, lb, nc - 1, la))
									        d1 = d1 / 2.0D0
									END IF

									! Application of the Wigner-Eckart theorem and multiplication by the mass factor
									! to get the proper energy
									d1 = d1 / (m(1) * DBLE(lla * llb))

									! Filling in the required pointers
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v11, la, na, nc, lb, nb, nd, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v11, lb, nb, nd, la, na, nc, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v22, la, na, nc, lb, nb, nd, cuad * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v22, lb, nb, nd, la, na, nc, cuad * d1)

									! Writing the results on disk
									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, "(6I3,2E24.16)", IOSTAT=file_error) &
										la, na, nc, lb, nb, nd, 0.5 * d1, cuad * d1
									END IF

									! PAIRING: < na la || NABLA || nb lb >< nc la || NABLA || nd lb >
									d1 = SymKineticEnergy2Body_nablaHO(na - 1, la, nb - 1, lb) * &
										SymKineticEnergy2Body_nablaHO(nc - 1, la, nd - 1, lb)

									! Anti-symmetrization of the matrix elements by permutation:
									!          < nb lb || NABLA || na la >< nd lb || NABLA || nc la >
									IF (nb .NE. nd) THEN
										d1 = d1 + (SymKineticEnergy2Body_nablaHO(na - 1, la, nd - 1, lb) * &
												SymKineticEnergy2Body_nablaHO(nc - 1, la, nb - 1, lb))
										d1 = d1 / 2.0D0
									END IF

									! Application of the Wigner-Eckart theorem and multiplication by the mass factor
									! to get the proper energy
									d1 = - d1 / (m(1) * DBLE(lla * llb))

									! Filling in the required pointers
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v11_pair, la, na, nc, lb, nb, nd, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v11_pair, lb, nb, nd, la, na, nc, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v22_pair, la, na, nc, lb, nb, nd, cuad * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v22_pair, lb, nb, nd, la, na, nc, cuad * d1)

								END DO
							END DO

						END DO

					END DO
				END DO

			END DO

		CASE (2)

			DO la = 0, Lmax

				lla = (2 * la) + 1 ! Degeneracy of a shell with l_a

							namax = MIN(Nmax, NmaxOfL(la))
				IF (CompHO .EQ. 1) 	namax = ((N_0 - la) / 2) + 1

				DO na = 1, namax
					DO nc = 1, na

						DO lb = 0, la

							llb = (2 * lb) + 1 ! Degeneracy of a shell with l_b

							cuad = CUAD2(la, lb, 1) ! Defined in Module math.f90

										nbmax = MIN(Nmax, NmaxOfL(lb))
							IF (CompHO .EQ. 1) 	nbmax = ((N_0 - lb) / 2) + 1

							DO nb = 1, nbmax
								DO nd = 1, nb

									! Product: < na la || NABLA || nd la >< nb lb || NABLA || nc lb >
									d1 = SymKineticEnergy2Body_nabla(na, la, nd, lb) * &
										SymKineticEnergy2Body_nabla(nb, lb, nc, la)

									! Anti-symmetrization of the matrix elements by permutation:
									!          < na la || NABLA || nb la >< nd lb || NABLA || nc lb >
									IF (nb .NE. nd) THEN
										d1 = d1 + (SymKineticEnergy2Body_nabla(na, la, nb, lb) * &
												SymKineticEnergy2Body_nabla(nd, lb, nc, la))
										d1 = d1 / 2.0D0
									END IF

									! Application of the Wigner-Eckart theorem and multiplication by the mass factor
									! to get the proper energy
									d1 = d1 / (m(1) * DBLE(lla * llb))

									! Filling in the required pointers
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v11, la, na, nc, lb, nb, nd, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v11, lb, nb, nd, la, na, nc, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v22, la, na, nc, lb, nb, nd, cuad * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v22, lb, nb, nd, la, na, nc, cuad * d1)

									! Writing the results on disk
									IF (file_error .EQ. 0) THEN
										WRITE (file_desc, "(6I3,2E24.16)", IOSTAT=file_error) &
												la, na, nc, lb, nb, nd, 0.5 * d1, cuad * d1
									END IF

									! PAIRING: < na la || NABLA || nb lb >< nc la || NABLA || nd lb >
									d1 = SymKineticEnergy2Body_nabla(na, la, nb, lb) * &
										SymKineticEnergy2Body_nabla(nc, la, nd, lb)

									! Anti-symmetrization of the matrix elements by permutation:
									!          < nb lb || NABLA || na la >< nd lb || NABLA || nc la >
									IF (nb .NE. nd) THEN
										d1 = d1 + (SymKineticEnergy2Body_nabla(na, la, nd, lb) * &
												SymKineticEnergy2Body_nabla(nc, la, nb, lb))
										d1 = d1 / 2.0D0
									END IF

									! Application of the Wigner-Eckart theorem and multiplication by the mass factor
									! to get the proper energy
									d1 = - d1 / (m(1) * DBLE(lla * llb))

									! Filling in the required pointers
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v11_pair, la, na, nc, lb, nb, nd, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v11_pair, lb, nb, nd, la, na, nc, 0.5D0 * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v22_pair, la, na, nc, lb, nb, nd, cuad * d1)
									CALL SymD3Tensor_SymD3Tensor_assign(vEkCMpp%v22_pair, lb, nb, nd, la, na, nc, cuad * d1)

								END DO
							END DO

						END DO

					END DO
				END DO

			END DO


		END SELECT

		CLOSE(file_desc)
		RETURN

	END SUBROUTINE SymKineticEnergy2Body_calculate

        !---------------------------------------------------------------------------------!
	!		                       		  			          !
	! Function giving the reduced matrix elements < n_a, l_a || NABLA || n_b, l_b >   !
	! for the special case of the harmonic oscillator basis                           !
	!		                       		  			          !
	! BEWARE: Here the HO length b is implicitely assumed equal to 1, otherwise,	  !
	!         there should be a factor 1/b						  !
	!		                       		  			          !
	! Refs: Appendix C, Page 115, Appendix D, Page 127, D24				  !
	!		                       		  			          !
        !---------------------------------------------------------------------------------!

	FUNCTION SymKineticEnergy2Body_nablaHO(na, la, nb, lb)
		DOUBLE PRECISION :: SymKineticEnergy2Body_nablaHO
		INTEGER, INTENT(IN) :: na, la, nb, lb

		IF (la .EQ. (lb + 1)) THEN
			IF (nb .EQ. (na + 1)) THEN
				SymKineticEnergy2Body_nablaHO = - sq(la) * sq(nb)/b_0
			ELSE IF (nb .EQ. na) THEN
				SymKineticEnergy2Body_nablaHO = - sq(la) * sq2(na + la)/b_0
			ELSE
				SymKineticEnergy2Body_nablaHO = 0.0D0
			END IF
		ELSE IF (lb .EQ. (la + 1)) THEN
			IF (na .EQ. (nb + 1)) THEN
				SymKineticEnergy2Body_nablaHO = - sq(lb) * sq(na)/b_0
			ELSE IF (na .EQ. nb) THEN
				SymKineticEnergy2Body_nablaHO = - sq(lb) * sq2(nb + lb)/b_0
			ELSE
				SymKineticEnergy2Body_nablaHO = 0.0D0
			END IF
		ELSE
			SymKineticEnergy2Body_nablaHO = 0.0D0
		END IF

		RETURN
	END FUNCTION SymKineticEnergy2Body_nablaHO

	!
	!  Calculating the kinetic energy fields by contracting the
	!  matrix elements with the density matrix
	!
	SUBROUTINE SymKineticEnergy2Body_get_Gamma(HF_out, vEkCMph, HF_in, ta, tb, ProjectionOn)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymKineticEnergy2Body), INTENT(IN) :: vEkCMph
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: HF_in
		INTEGER, INTENT(IN) :: ta, tb, ProjectionOn

		IF (ProjectionOn .EQ. 1) THEN

			IF (ta .EQ. tb) THEN
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(tb), vEkCMph%v11, HF_in%p(ta))
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(tb), vEkCMph%v22, HF_in%a(ta))
			ELSE
				HF_out%p(tb) = 0.0D0
				HF_out%a(tb) = 0.0D0
			END IF

		ELSE
			CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(tb), vEkCMph%v11, HF_in%p(ta))
			CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(tb), vEkCMph%v22, HF_in%a(ta))
		END IF

		HF_out%GaugeAngle(tb) = HF_in%GaugeAngle(ta)

		RETURN
	END SUBROUTINE SymKineticEnergy2Body_get_Gamma

	SUBROUTINE SymKineticEnergy2Body_get_Delta(HF_out, vEkCMpp, P_in, ta, tb, ProjectionOn)
		TYPE (SymHartreeFockFieldProj), INTENT(INOUT) :: HF_out
		TYPE (SymEk2pp), INTENT(IN) :: vEkCMpp
		TYPE (SymHartreeFockFieldProj), INTENT(IN) :: P_in
		INTEGER, INTENT(IN) :: ta, tb, ProjectionOn

		IF (ProjectionOn .EQ. 1) THEN

			IF (ta .EQ. tb) THEN
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(tb), vEkCMpp%v11_pair, P_in%p(ta))
				CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(tb), vEkCMpp%v22_pair, P_in%a(ta))
			ELSE
				HF_out%p(tb) = 0.0D0
				HF_out%a(tb) = 0.0D0
			END IF

		ELSE
			CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(tb), vEkCMpp%v11_pair, P_in%p(ta))
			CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(tb), vEkCMpp%v22_pair, P_in%a(ta))
		END IF

		HF_out%GaugeAngle(tb) = P_in%GaugeAngle(ta)

		RETURN
	END SUBROUTINE SymKineticEnergy2Body_get_Delta

	SUBROUTINE SymKineticEnergy2Body_del(vEkCMph, vEkCMpp)
		TYPE (SymKineticEnergy2Body), INTENT(INOUT) :: vEkCMph
		TYPE (SymEk2pp), INTENT(INOUT) :: vEkCMpp

		CALL SymD3Tensor_SymD3Tensor_del(vEkCMph%v11)
		CALL SymD3Tensor_SymD3Tensor_del(vEkCMph%v22)

		CALL SymD3Tensor_SymD3Tensor_del(vEkCMpp%v11_pair)
		CALL SymD3Tensor_SymD3Tensor_del(vEkCMpp%v22_pair)

		RETURN
	END SUBROUTINE SymKineticEnergy2Body_del

END MODULE symke2b
