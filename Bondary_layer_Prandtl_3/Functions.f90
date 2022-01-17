!************************************************************************************************      
      SUBROUTINE Solve_Tridiagonal_Matrix(NJ, A, B, C, D, y)
		IMPLICIT NONE
		
		integer :: J, NJ
		real(8), dimension(NJ) :: y
		real(8), dimension(NJ) :: A, B, C, D
		real(8), allocatable :: alpha(:), betta(:)
		real(8) :: temp
		allocate (alpha(NJ), betta(NJ))

        alpha(1) = -(C(1)/B(1))
        betta(1) = D(1)/B(1)
        do J=2, NJ
            temp = B(J) + A(J) * alpha(J-1)
            alpha(J) = -C(J) / temp
            betta(J) = (D(J) - A(J) * betta(J-1)) / temp
        enddo

        y(NJ) = betta(NJ)
        do J=NJ-1, 1, -1
            y(J) = alpha(J) * y(J+1) + betta(J)
        enddo
        deallocate (alpha,betta)

      END SUBROUTINE
	  
	  SUBROUTINE Find_tw(NI, NJ, x, y, tw, U, dy, MU)
		IMPLICIT NONE
		
		integer NI, NJ, i, j
		real(8), dimension(NI) :: tw
		
		real(8), dimension(NI, NJ) :: x, y, U
		
		real(8) dy, MU
		
		do i = 1, NI
			tw(i) = -MU * (3 * U(i, 1) - 4 * U(i, 2) + U(i, 3)) / (2 * dy)
		enddo
	  END SUBROUTINE