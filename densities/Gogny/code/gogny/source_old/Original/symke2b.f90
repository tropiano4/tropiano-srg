MODULE symke2b

	USE input
	USE global
	USE symd3t
	USE symfield

	IMPLICIT NONE

	TYPE SymKineticEnergy2Body
		TYPE (SymD3Tensor_SymD3Tensor) v11, v22
	END TYPE

!	INTERFACE OPERATOR(*)
!		MODULE PROCEDURE SymKineticEnergy2Body_product
!	END INTERFACE

	PRIVATE SymKineticEnergy2Body_nabla

CONTAINS

	SUBROUTINE SymKineticEnergy2Body_new(vEkCMph)
		TYPE (SymKineticEnergy2Body), INTENT(INOUT) :: vEkCMph

		INTEGER la, lla, na, namax, nc, lb, llb, nb, nbmax, nd
		DOUBLE PRECISION cuad, d1
		CHARACTER(64) filename

		INTEGER, PARAMETER :: file_desc = 6
		INTEGER file_error

		CALL SymD3Tensor_SymD3Tensor_new(vEkCMph%v11)
		CALL SymD3Tensor_SymD3Tensor_new(vEkCMph%v22)

		IF (N_0 < 10) THEN
			WRITE(filename, "(A,I1,A)") "data/vEkCM", N_0, "ph.txt"
		ELSE
			WRITE(filename, "(A,I2,A)") "data/vEkCM", N_0, "ph.txt"
		END IF

		IF (N_0 .LE. 16) THEN
			OPEN(file_desc, FILE=filename, ACTION="WRITE", IOSTAT=file_error)
			IF (file_error .NE. 0) THEN
				PRINT *, "*** Attention: Impossible to write the results in ", filename
			END IF
		END IF

		DO la = 0, N_0
			lla = (2 * la) + 1
			namax = ((N_0 - la) / 2) + 1
			DO na = 1, namax
!TODO			DO nc = 1, namax
! Debe haber sido un error del autor, puesto que estamos trabajando con
! matrices simétricas. Nota: el operador () corrige el error
				DO nc = 1, na
					DO lb = 0, la
						llb = (2 * lb) + 1
						nbmax = ((N_0 - lb) / 2) + 1
						cuad = CUAD2(la, lb, 1)
						DO nb = 1, nbmax
							DO nd = 1, nb
								d1 = SymKineticEnergy2Body_nabla(na - 1, la, nd - 1, lb) * &
								     SymKineticEnergy2Body_nabla(nb - 1, lb, nc - 1, la)
								IF (nb .NE. nd) THEN
									d1 = d1 + (SymKineticEnergy2Body_nabla(na - 1, la, nb - 1, lb) * &
									           SymKineticEnergy2Body_nabla(nd - 1, lb, nc - 1, la))
									d1 = d1 / 2.0
								END IF
								d1 = d1 / (m(1) * DBLE(lla * llb))

								CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v11, la, na, nc, lb, nb, nd, 0.5 * d1)
								CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v11, lb, nb, nd, la, na, nc, 0.5 * d1)
								CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v22, la, na, nc, lb, nb, nd, cuad * d1)
								CALL SymD3Tensor_SymD3Tensor_assign(vEkCMph%v22, lb, nb, nd, la, na, nc, cuad * d1)
					  IF ((N_0 .LE. 16) .AND. (file_error .EQ. 0)) THEN
					       WRITE (file_desc, "(I3,I3,I3,I3,I3,I3,E,E)", IOSTAT=file_error) &
								   la, na, nc, lb, nb, nd, 0.5 * d1, cuad * d1
					  END IF
								
							END DO
						END DO
					END DO
				END DO
			END DO
		END DO
		RETURN
	END SUBROUTINE SymKineticEnergy2Body_new

	FUNCTION SymKineticEnergy2Body_nabla(na, la, nb, lb)
		DOUBLE PRECISION SymKineticEnergy2Body_nabla
		INTEGER, INTENT(IN) :: na, la, nb, lb

		IF (la .EQ. (lb + 1)) THEN
			IF (nb .EQ. (na + 1)) THEN
				SymKineticEnergy2Body_nabla = - sq(la) * sq(nb)
			ELSE IF (nb .EQ. na) THEN
				SymKineticEnergy2Body_nabla = - sq(la) * sq2(na + la)
			ELSE
				SymKineticEnergy2Body_nabla = 0.0
			END IF
		ELSE IF (lb .EQ. (la + 1)) THEN
			IF (na .EQ. (nb + 1)) THEN
				SymKineticEnergy2Body_nabla = - sq(lb) * sq(na)
			ELSE IF (na .EQ. nb) THEN
				SymKineticEnergy2Body_nabla = - sq(lb) * sq2(nb + lb)
			ELSE
				SymKineticEnergy2Body_nabla = 0.0
			END IF
		ELSE
			SymKineticEnergy2Body_nabla = 0.0
		END IF
		RETURN
	END FUNCTION SymKineticEnergy2Body_nabla

	SUBROUTINE SymKineticEnergy2Body_product(HF_out, vEkCMph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymKineticEnergy2Body), INTENT(IN) :: vEkCMph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymKineticEnergy2Body_get_Gamma(HF_out, vEkCMph, HF_in)
		RETURN
	END SUBROUTINE SymKineticEnergy2Body_product

	SUBROUTINE SymKineticEnergy2Body_get_Gamma(HF_out, vEkCMph, HF_in)
		TYPE (SymHartreeFockField), INTENT(INOUT) :: HF_out
		TYPE (SymKineticEnergy2Body), INTENT(IN) :: vEkCMph
		TYPE (SymHartreeFockField), INTENT(IN) :: HF_in

		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(0), vEkCMph%v11, HF_in%p(0))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%p(1), vEkCMph%v11, HF_in%p(1))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(0), vEkCMph%v22, HF_in%a(0))
		CALL SymD3Tensor_SymD3Tensor_product(HF_out%a(1), vEkCMph%v22, HF_in%a(1))
		RETURN
	END SUBROUTINE SymKineticEnergy2Body_get_Gamma

	SUBROUTINE SymKineticEnergy2Body_del(vEkCMph)
		TYPE (SymKineticEnergy2Body), INTENT(INOUT) :: vEkCMph

		CALL SymD3Tensor_SymD3Tensor_del(vEkCMph%v11)
		CALL SymD3Tensor_SymD3Tensor_del(vEkCMph%v22)
		RETURN
	END SUBROUTINE SymKineticEnergy2Body_del

END MODULE symke2b
