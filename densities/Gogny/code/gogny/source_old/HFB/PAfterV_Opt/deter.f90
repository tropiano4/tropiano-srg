	
 MODULE deter
	
 CONTAINS
	
	!-------------------------------------------------------!
	!							!
	!     Nom :						!
	!       DETERM						!
	!							!
	!     Fonctionalite :					!
	!       Procedure de calcul de la valeur absolue du	!
	!       determinant					!
	!							!
	!     Entree(s) :					!
	!       A : matrice carree				!
	!       m : dimension de A				!
	!							!
	!     Sortie(s) :					!
	!       -> : valeur absolue du determinant		!
	!							!
	!-------------------------------------------------------!
	
	COMPLEX FUNCTION  DETERM (A, m)
	COMPLEX, DIMENSION(:,:), INTENT(IN) :: A
	INTEGER, INTENT(IN) :: m
	
	COMPLEX, DIMENSION(:,:), ALLOCATABLE :: B, C
	COMPLEX :: det
	INTEGER :: i, p, LIGNE
	
	i = 1
	
	ALLOCATE(B(m,m))
	ALLOCATE(C(m,m))

	B = A
	C = A
	
	DO WHILE (i .LT. m)
	
		p = LIGNE (B, m, i)
	
		IF (p .EQ. 0) THEN
			DETERM = CMPLX(0,0)
			RETURN
		END IF
	
		IF (i .NE. p) THEN
			CALL ECHANG (B, m, m, i, p, C)
			B = C
		END IF
	
		CALL PIVOTE (p, B, m, 1)
	
		i = i + 1
	 
	END DO 
	
	det = B(1, 1)
      
	DO i = 2, m
		det = det*B(i, i)
	END DO
 
	DEALLOCATE(B)
	DEALLOCATE(C)
	      
	DETERM = det
	
	RETURN
	END FUNCTION DETERM

	!-------------------------------------------------------!
	!							!
	!     Nom :						!
	!       LIGNE						!
	!							!
	!     Fonctionalite :					!
	!       Procedure de recherche de la ligne de pivot a	!
	!       l'etape n					!
	!							!
	!     Entree(s) :					!
	!       A : matrice carree				!
	!       m : dimension de A				!
	!       n : colonne de recherche			!
	!							!
	!     Sortie(s) :					!
	!       -> : ligne de pivot (0 si celle-ci n'existe pas)!
	!							!
	!-------------------------------------------------------!
	INTEGER FUNCTION LIGNE (A, m, n)
	COMPLEX, DIMENSION(:,:), INTENT(IN) :: A
	INTEGER, INTENT(IN) :: m, n
	
	INTEGER :: p
	
	p = n
 
 	DO WHILE ((p .LE. m) .AND. (ABS(A(p, n)) .LT. 1.E-10))
		p = p + 1
	END DO
      
	IF (p .GT. m) p = 0
	
	LIGNE = p
	
	RETURN
	END FUNCTION LIGNE


	!-------------------------------------------------------!
	!							!
	!     Nom :						!
	!      	 ECHANG						!
	!							!
	!     Fonctionalite :					!
	!       Procedure de permutation de deux ligne d'une 	!
	!       matrice ou d'un vecteur				!
	!							!
	!     Entree(s) :					!
	!       A : matrice / vecteur				!
	!       m, n : dimension				!
	!       k, l: lignes a permuter				!
	!							!
	!     Sortie(s) :					!
	!       B : matrice / vecteur				!
	!							!
	!-------------------------------------------------------!
	
	SUBROUTINE ECHANG (A, m, n, k, l, B)
	COMPLEX, DIMENSION(:,:), INTENT(OUT) :: B
	COMPLEX, DIMENSION(:,:), INTENT(IN) :: A
	INTEGER, INTENT(IN) :: m, n, k, l
	
	COMPLEX :: t
	INTEGER :: j
	
	B = A
	
	DO j = 1, n
		t = A(k, j)
		B(k, j) = A(l, j)
		B(l, j) = t
	END DO
 
	RETURN
	END SUBROUTINE ECHANG


	!-------------------------------------------------------!
	!							!
	!     Nom :						!
	!      	 PIVOT						!
	!							!
	!     Fonctionalite :					!
	!        Procedure de permettant d'effectuer le pivot a !
	!        l'etape i					!
	!							!
	!     Entree(s) :					!
	!       p : etape du pivot				!
	!       A : matrice carree				!
	!       m : dimension de A				!
	!       B: matrice / vecteur				!
	!       o : nombre de colonnes de B			!
	!							!
	!     Sortie(s) :					!
	!       A : matrice					!
	!       B : matrice / vecteur				!
	!							!
	!-------------------------------------------------------!
	
	SUBROUTINE PIVOTE (p, A, m, o)
	COMPLEX, DIMENSION(:,:), INTENT(OUT) :: A
	INTEGER, INTENT(IN) :: m, p, o
	
	INTEGER :: i, j
	COMPLEX :: piv, t
	
	piv = A(p, p)
	
	DO i = p+1, m
		
		t = A(i, p)
		A(i, p) = CMPLX(0,0)
		
		! Boucle de pivot de A
		DO j = p+1,m
			A(i, j) = A(i, j) - A(p, j)*t/piv
		END DO

	END DO
	
	RETURN
	END SUBROUTINE PIVOTE
	
 END MODULE deter
