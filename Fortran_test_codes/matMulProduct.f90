program matMulProduct

    integer, dimension(3,3) :: a, b, c
    integer :: i, j

    do i = 1, 3
        do j = 1, 3
            a(i, j) = i+j
        end do
    end do

    print *, 'Matrix Multiplication: A Matrix'

    do i = 1, 3
        do j = 1, 3
            print *, a(i, j)
        end do
    end do

    do i = 1, 3
        do j = 1, 3
            b(i, j) = i*j
        end do
    end do

    print *, 'Matrix Multiplication: B Matrix'

    do i = 1, 3
        do j = 1, 3
            print *, b(i, j)
        end do
    end do

    c = matmul(a, b)
    print *, 'Matrix Multiplication: Result Matrix'

    do i = 1, 3
        do j = 1, 3
            print *, c(i, j)
        end do
    end do

end program matMulProduct
