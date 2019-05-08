module magnus_wegner_functions
implicit none

    integer, dimension(2,2), private :: h0_matrix
    !real(DP), private :: delta = 1.0
    !real(DP), private :: g = 0.5
    integer, private :: n
    ! More stuff goes here

    ! Use test Hamiltonian from Hergert_2016iju
    !h0_matrix(1) = (/ 2.0*delta-g, -0.5*g, -0.5*g, -0.5*g, -0.5*g, 0.0 /)
    !h0_matrix(2) = (/ -0.5*g, 4.0*delta-g, -0.5*g, -0.5*g, 0.0, -0.5*g /)
    !h0_matrix(3) = (/ -0.5*g, -0.5*g, 6.0*delta-g, 0.0, -0.5*g, -0.5*g /)
    !h0_matrix(4) = (/ -0.5*g, -0.5*g, 0.0, 6.0*delta-g, -0.5*g, -0.5*g /)
    !h0_matrix(5) = (/ -0.5*g, 0.0, -0.5*g, -0.5*g, 8.0*delta-g, -0.5*g /)
    !h0_matrix(6) = (/ 0.0, -0.5*g, -0.5*g, -0.5*g, -0.5*g, 10.0*delta-g /)

    public :: commutator
    public :: diag

contains

    subroutine initialize_hamiltonian()

        integer :: i, j

        ! Do 2 x 2 test Hamiltonian

        !h0_matrix(1) = (/ 1, 1 /)
        !h0_matrix(2) = (/ 1, -1 /)

        do i = 1, n
            do j = 1, n
                if (i == n .and. j == n) then
                    h0_matrix(i, j) = -1
                else
                    h0_matrix(i, j) = 1
                end if
            end do
        end do

        ! Dimension of Hamiltonian
        n = int( sqrt( real( size(h0_matrix) ) ) )

    end subroutine initialize_hamiltonian

    function commutator(matrix_1, matrix_2)result(matrix_3)

        integer, dimension(n,n), intent(in) :: matrix_1, matrix_2
        integer, dimension(n,n) :: matrix_3

        matrix_3 = matmul(matrix_1, matrix_2) - matmul(matrix_2, matrix_1)

    end function commutator

    function diag(matrix_1)result(matrix_2)

        integer :: i, j
        integer, dimension(n,n), intent(in) :: matrix_1
        integer, dimension(n,n) :: matrix_2

        do i = 1, n
            do j = 1, n
                if (i == j) then
                    matrix_2(i, j) = matrix_1(i, j)
                else
                    matrix_2(i, j) = 0
                end if
            end do
        end do

    end function diag


end module magnus_wegner_functions

program magnus_wegner
use magnus_wegner_functions
implicit none

    integer, dimension(n,n) :: hd_matrix
    integer, dimension(n,n) :: eta

    hd_matrix = diag(h0_matrix)
    eta = commutator(hd_matrix, h0_matrix)

    print *, 'Eta: Result Matrix'

    do i = 1, n
        do j = 1, n
            print *, eta(i, j)
        end do
    end do

end program magnus_wegner
