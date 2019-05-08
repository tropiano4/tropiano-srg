module magnus_wegner
implicit none


	! Dimension of the Hamiltonian
	integer, parameter :: n = 6
	! Parameters for test pairing Hamiltonian
	double precision, parameter :: delta = 0.1, g = 0.05
	! Number of terms to include in the Magnus sum
	integer, parameter :: k_magnus = 7
	! Declare initial Hamiltonian
	double precision, dimension(n, n) :: hamiltonian_0
	! Declare factorial and Bernoulli number arrays
	double precision, dimension(k_magnus) :: factorial_array, bernoulli_array, magnus_factors
	! Step-size in solving omega differential equation
	double precision, parameter :: ds = 1.0e-5

	! put public functions or subroutines here


contains


	subroutine initialize()
		! This subroutine initializes the bare Hamiltonian and arrays for storing
		! factorials and Bernoulli numbers.

		integer :: i, j, k

		! Use 6 x 6 Hamiltonian as test case from Hergert_2016iju
		print *, 'Initial hamiltonian'
		do i = 1, n
			do j = 1, n
				if (i == j) then
					if (i > 3) then
						hamiltonian_0(i, j) = 2.0*dble(i-1)*delta - g
					else
						hamiltonian_0(i, j) = 2.0*dble(i)*delta - g
					end if
				else
					if ( (i+j) == 7 ) then
						hamiltonian_0(i, j) = 0.0
					else
						hamiltonian_0(i, j) = -0.5*g
					end if
				end if
				print *, hamiltonian_0(i, j)
			end do
		end do

		! Factorial and Bernoulli numbers up through k_magnus
		factorial_array = (/ 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0 /)
		bernoulli_array = (/ 1.0, -1.0/2.0, 1.0/6.0, 0.0, -1.0/30.0, 0.0, 1.0/42.0 /)

		! Factors for Magnus sum
		do k = 1, k_magnus
			magnus_factors(k) = bernoulli_array(k) / factorial_array(k)
		end do

	end subroutine initialize


	function commutator(matrix_1, matrix_2)result(matrix_3)
		! Returns the commutator of two matrices.

		double precision, dimension(n, n), intent(in) :: matrix_1, matrix_2
		double precision, dimension(n, n) :: matrix_3

		matrix_3 = matmul(matrix_1, matrix_2) - matmul(matrix_2, matrix_1)

	end function commutator


	function diag(matrix_1)result(matrix_2)
		! Returns a matrix of the diagonal elements of the given matrix.

		integer :: i, j
		double precision, dimension(n, n), intent(in) :: matrix_1
		double precision, dimension(n, n) :: matrix_2

		do i = 1, n
			do j = 1, n
				! Keep diagonal elements
				if (i == j) then
					matrix_2(i, j) = matrix_1(i, j)
				! Set off-diagonal elements to zero
				else
					matrix_2(i, j) = 0.0
				end if
			end do
		end do

	end function diag


	function bch_formula(operator_0, omega, k_terms)result(operator_s)
		! Returns an evolved operator given an initial operator and evolved
		! omega. Must specify truncation.

		! Initial operator matrix
		double precision, dimension(n, n), intent(in) :: operator_0
	 	! Evolving omega matrix
		double precision, dimension(n, n), intent(in) :: omega
		! Number of terms to include in the BCH sum
		integer, intent(in) :: k_terms
		integer :: k
		! ad_omega^k ( operator ) - see Magnus paper for details
		double precision, dimension(n, n) :: ad
		! Evolved operator matrix
		double precision, dimension(n, n) :: operator_s

		! Initial nested commutator ad(operator)
		ad = operator_0

		! Zeroth term of the evolved operator
		operator_s = ad/factorial_array(1)

		! Sum through 1 to k_terms (note k = 1 corresponds to the 0th term and so forth)
		do k = 2, k_terms

			ad = commutator(omega, ad)
			operator_s = operator_s + ad/factorial_array(k)

		end do

	end function bch_formula


	function derivs(omega)result(domega)
		! Returns the RHS of the Magnus Omega equation where omega is the
        ! solution matrix.'''

		! Evolving omega matrix
		double precision, dimension(n, n), intent(in) :: omega
		integer :: k
		! Evolving Hamiltonian matrix
		double precision, dimension(n, n) :: hs_matrix
		! G matrix for Wegner SRG generator
		double precision, dimension(n, n) :: hd_matrix
		! eta SRG generator
		double precision, dimension(n, n) :: eta
		! ad_omega^k ( eta ) - see Magnus paper for details
		double precision, dimension(n, n) :: ad
		! Derivative of the evolving omega matrix
		double precision, dimension(n, n) :: domega

		! Obtain evolved Hamiltonian with BCH formula summing through 25 terms
		hs_matrix = bch_formula(hamiltonian_0, omega, k_magnus)

		! Wegner SRG generator, eta = [G,H] where G = H_D
		hd_matrix = diag(hs_matrix)

		! SRG generator [G, H(s)]
		eta = commutator(hd_matrix, hs_matrix)

		! Initial nested commutator ad(eta)
		ad = eta

		domega = ad*magnus_factors(1)

		! Sum through 1 to k_magnus (note k = 1 corresponds to the 0th term and so forth)
		do k = 2, k_magnus

			ad = commutator(omega, ad)
			domega = domega + magnus_factors(2)*ad

		end do

	end function derivs


	function euler_method(omega_0, s_init, s_max)result(omega_s)
		! Use first-order Euler method with fixed step-size ds to solve Magnus
		! Omega equation.

		! Initial omega matrix
		double precision, dimension(n, n), intent(in) :: omega_0
		! Initial value of s
		double precision, intent(in) :: s_init
		! Maximum value of s
		double precision, intent(in) :: s_max
		! Flowing value of s in the range s_init to s_max
		double precision :: s
		! Evolved omega matrix
		double precision, dimension(n, n) :: omega_s

		omega_s = omega_0

		! Step through s values until fully evolved
		s = s_init

		do while ( s <= s_max )

			omega_s = omega_s + derivs(omega_s)
			s = s + ds

		end do

	end function euler_method


	function evolve_hamiltonian(lambda)result(hamiltonian_s)
		! Returns evolved Hamiltonian at the given value of lambda.

		! Lambda evolution value in units fm^-1
		double precision, intent(in) :: lambda
		integer :: i, j
		! Convert from lambda to s
		double precision :: s_init = 0.0
		double precision :: s_max
		! Initial and evolved omega matrices
		double precision, dimension(n, n) :: omega_init, omega_evolved
		! Evolved Hamiltonian matrix
		double precision, dimension(n, n) :: hamiltonian_s

		s_max = 1.0/lambda**4.0

		! Initial omega matrix is entirely zero
		do i = 1, n
			do j = 1, n
				omega_init(i, j) = 0.0
			end do
		end do

		! Evolve omega with first-order Euler method solver
		omega_evolved = euler_method(omega_init, s_init, s_max)

		! Evolve Hamiltonian using the BCH formula
		hamiltonian_s = bch_formula(hamiltonian_0, omega_evolved, k_magnus)

	end function evolve_hamiltonian


end module magnus_wegner


program test_magnus
use magnus_wegner
implicit none


	! Declare variables
	integer :: i, j
	! Lambda evolution value in units fm^-1
	double precision :: lamb = 3.0
	! Declare evolved Hamiltonian matrix
	double precision, dimension(n, n) :: hamiltonian_s

	! Initialize Hamiltonian and constants
	call initialize()

	hamiltonian_s = evolve_hamiltonian(lamb)

	print *, 'H(s): Result Matrix'

	do i = 1, n
		do j = 1, n
			print *, hamiltonian_s(i, j)
		end do
	end do


end program test_magnus
