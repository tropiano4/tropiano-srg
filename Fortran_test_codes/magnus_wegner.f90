module magnus_wegner_functions
implicit none

    !real(DP), dimension(6,6), private :: h0_matrix
    !real(DP), private :: delta = 1.0
    !real(DP), private :: g = 0.5
    integer, private :: n = 3
    ! More stuff goes here

    ! Use test Hamiltonian from Hergert_2016iju
    !h0_matrix(1) = (/ 2.0*delta-g, -0.5*g, -0.5*g, -0.5*g, -0.5*g, 0.0 /)
    !h0_matrix(2) = (/ -0.5*g, 4.0*delta-g, -0.5*g, -0.5*g, 0.0, -0.5*g /)
    !h0_matrix(3) = (/ -0.5*g, -0.5*g, 6.0*delta-g, 0.0, -0.5*g, -0.5*g /)
    !h0_matrix(4) = (/ -0.5*g, -0.5*g, 0.0, 6.0*delta-g, -0.5*g, -0.5*g /)
    !h0_matrix(5) = (/ -0.5*g, 0.0, -0.5*g, -0.5*g, 8.0*delta-g, -0.5*g /)
    !h0_matrix(6) = (/ 0.0, -0.5*g, -0.5*g, -0.5*g, -0.5*g, 10.0*delta-g /)

    ! Dimension of Hamiltonian
    !n = size(h0_matrix(1))

    public :: commutator

contains

    function commutator(matrix_1, matrix_2)result(matrix_3)

        !real(DP), dimension(n,n), private :: matrix_1, matrix_2, matrix_3
        integer, dimension(n,n), intent(in) :: matrix_1, matrix_2
        integer, dimension(n,n) :: matrix_3

        matrix_3 = matmul(matrix_1, matrix_2) - matmul(matrix_2, matrix_1)

    end function commutator

end module magnus_wegner_functions

program magnus_wegner
use magnus_wegner_functions
implicit none

    integer, dimension(3,3) :: a, b, c
    integer :: i, j

    do i = 1, 3
        do j = 1, 3
            a(i, j) = i+j
        end do
    end do

    do i = 1, 3
        do j = 1, 3
            b(i, j) = i*j
        end do
    end do

    c = commutator(a, b)
    print *, 'Commutator of A and B: Result Matrix'

    do i = 1, 3
        do j = 1, 3
            print *, c(i, j)
        end do
    end do

end program magnus_wegner
