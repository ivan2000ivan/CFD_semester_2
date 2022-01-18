
!************************************************************************************************      
      SUBROUTINE Prandtl(NI, NJ, X, Y, U, V, P, R, EPS, SMAX, dx, dy, NU, H, R0, P0)
		IMPLICIT NONE
		
		INTEGER NI,NJ, i, j, Su, Sp, SMAX
        real(8),DIMENSION(NI,NJ):: X,Y
        real(8),DIMENSION(NI, NJ):: U,V,P, R
		real(8) EPS, dx, dy, NU, e_U, e_U1, e_V, e_V1, H, R0, P0, Const, e_P, e_P1, betta
		
		real(8),DIMENSION(NJ):: Uold, Vold, Pold, Rold

		real(8),ALLOCATABLE :: A(:), B(:), C(:), D(:), G(:)
		allocate(A(NJ), B(NJ), C(NJ), D(NJ), G(NI))

		G = 0
		do j = 2, NJ
			G(1) = G(1) + (R(1, j) * U(1, j) + R(1, j - 1) * U(1, j - 1)) * dy / 2.0
		enddo
		
		!G(1) = R0 * 90 * 0.1
		
		Const = P0 / R0**1.4
		betta = -G(1) / (R0 * H**2)
		
		A(1) = 0.
		B(1) = -1.
		C(1) = 1.
		D(1) = 0.!U(i, 2)
		
		A(NJ) = -1.
		B(NJ) = 1.
		C(NJ) = 0.
		D(NJ) = 0.!U(i, NJ - 1)
		
		
		do i = 2, NI 
			Sp = 0
			P(i, :) = P(i - 1, :)
			R(i, :) = (P(i, :) / Const)**(1.0/1.4)
			U(i, :) = U(i - 1, :)
			V(i, :) = V(i - 1, :)
			e_P = 1
			do while(e_P > 1.0E-6 )!.AND. Sp < SMAX)
				!P(i, :) = P(i - 1, :)
				Su = 0
				Pold(1:NJ) = P(i, 1:NJ)
				Rold(1:NJ) = R(i, 1:NJ)
				!U(i, :) = U(i - 1, :)
				!V(i, :) = V(i - 1, :)
				e_V = 1
				do while (max(e_U, e_V) > 1.0E-6 )!.AND. Su < SMAX)  !(Su < SMAX)
					!U(i, :) = U(i - 1, :)
					!V(i, :) = V(i - 1, :)
					Uold(1:NJ) = U(i, 1:NJ)
					Vold(1:NJ) = V(i, 1:NJ)
					do j = 2, NJ - 1
						A(j) = -((V(i, j - 1) / (2.0 * dy)) + (NU / (dy)**2)) * R(i, j) !+
						B(j) = ((U(i, j) / dx) + (2.0 * NU / (dy)**2)) * R(i, j)
						C(j) = ((V(i, j + 1) / (2.0 * dy)) - (NU / (dy)**2)) * R(i, j)
						D(j) = (R(i - 1, j) * (U(i - 1, j))**2) / dx - (P(i, j) - P(i - 1, j)) / dx
					enddo
					
					
					call Solve_Tridiagonal_Matrix(NJ, A, B, C, D, U(i, 1:NJ))
				
					V(i, 1) = 0.
					V(i, NJ) = 0.
					do j = 2, NJ - 1
						V(i, j) = V(i, j-1) * R(i, j) - dy * (U(i, j) * R(i, j) - U(i - 1, j) * R(i-1, j) + U(i, j - 1) * R(i, j) - U(i - 1, j - 1) * R(i-1, j)) / (2 * dx)
						V(i, j) = V(i, j) / R(i, j)
					enddo
					
					Su = Su + 1
					
					e_U1 = abs(U(i, 1))
					e_V1 = abs(V(i, 1))
					do j = 2, NJ
						e_U1 = max(e_U1, abs(U(i, j)))
						e_V1 = max(e_V1, abs(V(i, j)))
					enddo
					e_U = abs(U(i, 1) - Uold(1)) / e_U1
					e_V = abs(V(i, 1) - Vold(1)) / e_V1
					do j = 2, NJ
						e_U = max(e_U, abs(U(i, j) - Uold(j)) / e_U1)
						e_V = max(e_V, abs(V(i, j) - Vold(j)) / e_V1)
					enddo
					!print*, "SU = ", Su, e_U, e_V, e_U1, e_V1
				enddo
				!pause
				Sp = Sp + 1
				G(i) = 0
				do j = 2, NJ
					G(i) = G(i) + ((R(i, j) * U(i, j) + R(i, j - 1) * U(i, j - 1)) * dy / 2.0)
				enddo
				!betta = -G(i) / (R0 * H**2)
				do j = 1, NJ
					betta = -G(i) / (R(i, j) * H**2)
					P(i, j) = P(i, j) + betta * (G(i - 1) - G(i))
					R(i, j) = (P(i, j) / Const)**(1.0/1.4)
				enddo
				print*, i, P(i, 5), R(i, 5), G(i), betta * (G(i - 1) - G(i))
				!pause
				e_P1 = P(i, 1)
				do j = 2, NJ
					e_P1 = max(e_P1, P(i, j))
				enddo
				e_P = abs(P(i, 1) - Pold(1)) / e_P1
				do j = 2, NJ
					e_P = max(e_P, abs(P(i, j) - Pold(j)) / e_P1)
				enddo
				print*, i, "SP = ", Sp, e_P
				!pause
			enddo
		enddo
		
		deallocate(A, B, C, D, G)
      END SUBROUTINE
	  
