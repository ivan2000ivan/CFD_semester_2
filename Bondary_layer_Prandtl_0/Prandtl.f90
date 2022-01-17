
!************************************************************************************************      
      SUBROUTINE Prandtl(NI, NJ, X, Y, U, V, P, U0, EPS, SMAX, dx, dy, NU)
		IMPLICIT NONE
		
		INTEGER NI,NJ, i, j, s, SMAX
        REAL,DIMENSION(NI,NJ):: X,Y
        REAL,DIMENSION(NI,NJ):: U,V,P
		Real U0, EPS, dx, dy, NU, e, e1
		
		REAL,DIMENSION(NJ):: Uold

		REAL,ALLOCATABLE :: A(:), B(:), C(:), D(:)
		allocate(A(NJ), B(NJ), C(NJ), D(NJ))
		
		A(1) = 0
		B(1) = 1
		C(1) = 0
		D(1) = 0
		
		A(NJ) = 0
		B(NJ) = 1
		C(NJ) = 0
		D(NJ) = U0
		
		!s = 1
		e = 1.0e-6
		e1 = 1.0e-6
		
		!do while(s < SMAX)
			do i = 2, NI 
				s = 1
				Uold(1:NJ) = V(i, 1:NJ)
				do while(s < SMAX)
					Uold(1:NJ) = U(i, 1:NJ)
					do j = 2, NJ - 1
				
						A(j) = -((V(i, j - 1) / (2 * dy)) + (NU / (dy)**2))
						B(j) = (U(i, j) / dx) + (2 * NU / (dy)**2)
						C(j) = (V(i, j + 1) / (2 * dy)) - (NU / (dy)**2)
						D(j) = ((U(i - 1, j))**2) / dx
					enddo
					!print*, A
					!print*, "--------"
					!print*, B
					!print*, "--------"
					!print*, C
					!print*, "--------"
					!print*, D
					!print*, "--------"
				
					call Solve_Tridiagonal_Matrix(NJ, U(i, 1:NJ), C, B, A, -D)
				
					do j = 2, NJ
						V(i, j) = V(i, j-1) - dy * (U(i, j) - U(i - 1, j) + U(i, j - 1) - U(i - 1, j - 1)) / (2 * dx)
					enddo
					s = s + 1
					
					e = 1.0e-6
					e1 = 1.0e-6
					do j = 1, NJ
						e1 = max(e1, U(i, j))
					enddo
					do j = 1, NJ
						e = max(e, abs(U(i, j) - Uold(j)) / e1)
					enddo
					print*, e
				enddo
			enddo
			
		!enddo
		
		
		deallocate(A, B, C, D)
      END SUBROUTINE
	  
