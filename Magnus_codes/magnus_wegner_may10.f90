module magnus_wegner

	! Dimension of the Hamiltonian
	integer, parameter :: n = 240
	! Number of terms to include in the Magnus sum
	integer, parameter :: k_magnus = 7
	! Declare initial Hamiltonian
	double precision, dimension(n, n) :: hamiltonian_0
	! Declare factorial and Bernoulli number arrays
	double precision, dimension(k_magnus) :: factorial_array, bernoulli_array, magnus_factors
	! Step-size in solving omega differential equation
	double precision, parameter :: ds = 1.0e-5


contains


	subroutine initialize()
		! This subroutine arrays for storing factorials and Bernoulli numbers.

		integer :: k

		! Factorial and Bernoulli numbers up through k_magnus
		factorial_array = (/ 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0 /)
		bernoulli_array = (/ 1.0, -1.0/2.0, 1.0/6.0, 0.0, -1.0/30.0, 0.0, 1.0/42.0 /)

		! Factors for Magnus sum
		do k = 1, k_magnus
			magnus_factors(k) = bernoulli_array(k) / factorial_array(k)
		end do

	end subroutine initialize


	function commutator(matrix_1, matrix_2) result(matrix_3)
		! Returns the commutator of two matrices.

		double precision, dimension(240, 240) :: matrix_1, matrix_2
		double precision, dimension(240, 240) :: matrix_3

		matrix_3 = matmul(matrix_1, matrix_2) - matmul(matrix_2, matrix_1)

	end function commutator


	function diag(matrix_1) result(matrix_2)
		! Returns a matrix of the diagonal elements of the given matrix.

		integer :: i, j
		double precision, dimension(240, 240) :: matrix_1
		double precision, dimension(240, 240) :: matrix_2

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


	function bch_formula(operator_0, omega, k_terms) result(operator_s)
		! Returns an evolved operator given an initial operator and evolved
		! omega. Must specify truncation.

		! Initial operator matrix
		double precision, dimension(240, 240) :: operator_0
		! Evolving omega matrix
		double precision, dimension(240, 240) :: omega
		! Number of terms to include in the BCH sum
		integer :: k_terms
		integer :: k
		! ad_omega^k ( operator ) - see Magnus paper for details
		double precision, dimension(240, 240) :: ad
		! Evolved operator matrix
		double precision, dimension(240, 240) :: operator_s

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


	function derivs(omega) result(domega)
		! Returns the RHS of the Magnus Omega equation where omega is the
		! solution matrix.'''

		! Evolving omega matrix
		double precision, dimension(240, 240) :: omega
		integer :: k
		! Evolving Hamiltonian matrix
		double precision, dimension(240, 240) :: hs_matrix
		! G matrix for Wegner SRG generator
		double precision, dimension(240, 240) :: hd_matrix
		! eta SRG generator
		double precision, dimension(240, 240) :: eta
		! ad_omega^k ( eta ) - see Magnus paper for details
		double precision, dimension(240, 240) :: ad
		! Derivative of the evolving omega matrix
		double precision, dimension(240, 240) :: domega

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


	function euler_method(omega_0, s_init, s_max) result(omega_s)
		! Use first-order Euler method with fixed step-size ds to solve Magnus
		! Omega equation.

		! Initial omega matrix
		double precision, dimension(240, 240) :: omega_0
		! Initial value of s
		double precision :: s_init
		! Maximum value of s
		double precision :: s_max
		! Flowing value of s in the range s_init to s_max
		double precision :: s
		! Evolved omega matrix
		double precision, dimension(240, 240) :: omega_s

		omega_s = omega_0

		! Step through s values until fully evolved
		s = s_init

		do while ( s <= s_max )

			omega_s = omega_s + derivs(omega_s)
			s = s + ds

		end do

	end function euler_method


	function evolve_hamiltonian(h0_matrix, lambda) result(hamiltonian_s)
		! Returns evolved Hamiltonian at the given value of lambda.

		! h-bar ^ 2 / M [ MeV fm^2 ]
		double precision :: hbar_sq_over_M = 41.47
		! Bare Hamiltonian matrix
		double precision, dimension(240, 240) :: h0_matrix
		! Lambda evolution value in units fm^-1
		double precision :: lambda
		integer :: i, j
		! Convert from lambda to s
		double precision :: s_init = 0.0
		double precision :: s_max
		! Initial and evolved omega matrices
		double precision, dimension(240, 240) :: omega_init, omega_evolved
		! Evolved Hamiltonian matrix
		double precision, dimension(240, 240) :: hamiltonian_s

		! f2py intent(in) :: h0_matrix
		! f2py intent(in) :: lambda
		! f2py intent(out) :: hamiltonian_s

		! Convert Hamiltonian to fm^-2
		hamiltonian_0 = h0_matrix / hbar_sq_over_M

		! Initialize constants
		call initialize()

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
