 MODULE indexx

 CONTAINS

	SUBROUTINE indexx_real8(n, arrin, indx)
		INTEGER, INTENT(IN) :: n
		DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: arrin
		INTEGER, DIMENSION(*), INTENT(INOUT) :: indx

		INTEGER l, j, ir, indxt, i
		DOUBLE PRECISION q

		DO j = 1, n
			indx(j) = j
		END DO

		IF (n .EQ. 1) RETURN

		l = (n / 2) + 1
		ir = n
		DO
			IF (l .GT. 1) THEN
				l = l - 1
				indxt = indx(l)
				q = arrin(indxt)
			ELSE
				indxt = indx(ir)
				q = arrin(indxt)
				indx(ir) = indx(1)
				ir = ir - 1
				IF (ir .EQ. 1) THEN
					indx(1) = indxt
					RETURN
				END IF
			END IF
			i = l
			j = l * 2
			DO WHILE (j .LE. ir)
			
				IF (j .LT. ir) THEN
					IF (arrin(indx(j)) .LT. arrin(indx(j + 1))) THEN
						j = j + 1
					END IF
				END IF
				
				IF (q .LT. arrin(indx(j))) THEN
					indx(i) = indx(j)
					i = j
					j = j + i
				ELSE
					j = ir + 1
				END IF
			END DO
			indx(i) = indxt
		END DO
	END SUBROUTINE indexx_real8

 END MODULE indexx
