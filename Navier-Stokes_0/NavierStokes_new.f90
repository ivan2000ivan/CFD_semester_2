!************************************************************************************************      
       SUBROUTINE NavierStokes(x, y, U_c, V_c, P_c, dx, dy, U0, NI, NJ, NU, EPS, CFL, NITER)
       IMPLICIT NONE
	   
	   INTEGER NI, NJ, NITER, i, j, k
	   REAL(8),DIMENSION(0:NI,0:NJ):: Ug, Vg, Uv, Vv, Pg, Pv, U_c, V_c, P_c, Uold, Vold, Pold, Resp, Resx, Resy
	   REAL(8),DIMENSION(NI,NJ):: x, y
	   REAL(8) dx, dy, EPS, CFL, dt, A, Ur, Ul, Ut, Ub, Vt, Vb, Vr, Vl, Pl, Pr, Pt, Pb, NU, U0, Ru, Rv, Rp, MASS, IMPLX, IMPLY, Ue, Ve
	   
	   INTEGER IU, IV, IP
	   
	   IU = 1
	   IV = 2
	   IP = 3
	   
	   open(IU, file = "Ru.txt")
	   open(IV, file = "Rv.txt")
	   open(IP, file = "Rp.txt")
	   
	   dt = CFL * (dx + dy) / (2.0 * U0)
	   
	   A = 1 / (U0**2)
	   
	   Resp = 0
	   Resx = 0
	   Resy = 0

	   
	   do k = 1, NITER
		Uold = U_c
		Vold = V_c
		Pold = P_c
		Ru = 0
		Rv = 0
		Rp = 0
		
		do j = 1, NJ - 1
			do i = 0, NI - 1
				Ue = (U_c(i + 1, j) + U_c(i, j)) * 0.5
				if (Ue .GE. 0) then
					Ut = U_c(i, j)
					Pt = P_c(i + 1, j)
					Vt = V_c(i, j)
				else
					Ut = U_c(i + 1, j)
					Pt = P_c(i, j)
					Vt = V_c(i + 1, j)
				endif
				MASS = Ut
				IMPLX = Ue * Ut + Pt - NU * (U_c(i + 1, j) - U_c(i, j)) / dx
				IMPLY = Ue * Vt - NU * (V_c(i + 1, j) - V_c(i, j)) / dx
				
				Resp(i, j) = Resp(i, j) + MASS / dx
				Resp(i + 1, j) = Resp(i + 1, j) - MASS / dx
				
				Resx(i, j) = Resx(i, j) + IMPLX / dx
				Resx(i + 1, j) = Resx(i + 1, j) - IMPLX / dx
				
				Resy(i, j) = Resy(i, j) + IMPLY / dx
				Resy(i + 1, j) = Resy(i + 1, j) - IMPLY / dx
			enddo
		enddo
		
		do i = 1, NI - 1
			do j = 0, NJ - 1
				Ve = (V_c(i, j + 1) + V_c(i, j)) * 0.5
				if (Ve .GE. 0) then
					Ut = U_c(i, j)
					Pt = P_c(i, j + 1)
					Vt = V_c(i, j)
				else
					Ut = U_c(i, j + 1)
					Pt = P_c(i, j)
					Vt = V_c(i, j + 1)
				endif
				MASS = Vt
				IMPLX = Ve * Ut - NU * (U_c(i, j + 1) - U_c(i, j)) / dy
				IMPLY = Ve * Vt + Pt - NU * (V_c(i, j + 1) - V_c(i, j)) / dy
				
				Resp(i, j) = Resp(i, j) + MASS / dy
				Resp(i, j + 1) = Resp(i, j + 1) - MASS / dy
				
				Resx(i, j) = Resx(i, j) + IMPLX / dy
				Resx(i, j + 1) = Resx(i, j + 1) - IMPLX / dy
				
				Resy(i, j) = Resy(i, j) + IMPLY / dy
				Resy(i, j + 1) = Resy(i, j + 1) - IMPLY / dy
			enddo
		enddo
		
		Resp = -dt * U0**2 * Resp
		Resx = -dt * Resx
		Resy = -dt * Resy
		
	    P_c = P_c + Resp
		U_c = U_c + Resx
		V_c = V_c + Resy
		
		
		
		!ГУ 
		
		!Непроницаемая граница:
		U_c(0:NI, 0) = -U_c(0:NI, 1)
		V_c(0:NI, 0) = -V_c(0:NI, 1)
		P_c(0:NI, 0) = P_c(0:NI, 1)
		
		!Входная граница:
		U_c(0, 1:NJ) = U0
		V_c(0, 1:NJ) = 0.0
		P_c(0, 1:NJ-1) = P_c(1, 1:NJ-1)
		
		!Симметрия по х и у:
		
		!Верхнняя граница:
		do i = 1, NI - 1
			V_c(i, NJ) = V_c(i, NJ - 1)
			if(V_c(i, NJ - 1) > 0) then
				U_c(i, NJ) = U_c(i, NJ - 1)
				P_c(i, NJ) = 0.0
			else
				U_c(i, NJ) = U0
				P_c(i, NJ) = P_c(i, NJ-1)
			endif
		enddo
		
		!Выходная граница:
		U_c(NI, 0:NJ) = U_c(NI - 1, 0:NJ)
		V_c(NI, 0:NJ) = V_c(NI - 1, 0:NJ)
		P_c(NI, 0:NJ) = 0.0
		
		!Невязки:
		
		do i = 1, NI - 1
			do j = 1, NJ - 1
				Ru = max(abs(Resx(i, j)) / dt, Ru)
				Rv = max(abs(Resy(i, j)) / dt, Rv)
				Rp = max(abs(Resp(i, j)) / dt, Rp)
			enddo
		enddo
		
		write(IU, *) k, Ru
		write(IV, *) k, Rv
		write(IP, *) k, Rp
		
	   enddo
	   
	   close(IU)
	   close(IV)
	   close(IP)

       END  SUBROUTINE

!************************************************************************************************
                    