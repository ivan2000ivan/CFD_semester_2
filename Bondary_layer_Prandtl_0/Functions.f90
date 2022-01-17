!************************************************************************************************      
      SUBROUTINE Solve_Tridiagonal_Matrix(N, y, A, B, C, D)
		IMPLICIT NONE
		integer i, Im, N
		real(4), dimension(N) :: y, A, B, C, D, alfa, beta
        real E
        Im = size(A)
        !allocate (alfa(Im), beta(Im))
		
        alfa(1) = -A(1) / B(1)
        beta(1) = -D(1) / B(1)
        do i = 2, Im - 1
            E = B(i) + C(i) * alfa(i - 1)
            alfa(i) = -A(i) / E
            beta(i) = -(D(i) + C(i) * beta(i - 1)) / E
        enddo
        y(Im) = -(D(Im) + C(Im) * beta(Im - 1)) / (B(Im) + C(Im) * alfa(Im - 1))
        do i = Im - 1, 1, -1
            y(i) = alfa(i) * y(i + 1) + beta(i)
        enddo
        !deallocate(alfa, beta)
      END SUBROUTINE
	  
	  SUBROUTINE Find_tw(NI, NJ, x, y, tw, U, dy, MU)
		IMPLICIT NONE
		
		integer NI, NJ, i, j
		real(4), dimension(NI) :: tw
		
		real(4), dimension(NI, NJ) :: x, y, U
		
		real dy, MU
		
		do i = 1, NI
			tw(i) = -MU * (3 * U(i, 1) - 4 * U(i, 2) + U(i, 3)) / (2 * dy)
		enddo
	  END SUBROUTINE