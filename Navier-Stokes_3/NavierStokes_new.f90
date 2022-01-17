!************************************************************************************************      
       SUBROUTINE NavierStokes(x, y, U_c, V_c, P_c, R_c, dx, dy, U0, NI, NJ, MU, EPS, CFL, NITER, P0, R0, h_V, U1, U2, Y_c)
       IMPLICIT NONE
	   
	   INTEGER NI, NJ, NITER, i, j, k
	   REAL(8),DIMENSION(0:NI,0:NJ):: Ug, Vg, Uv, Vv, Pg, Pv, U_c, V_c, P_c, Uold, Vold, Pold, Resp, Resx, Resy, R_c, Rold, Y_c
	   REAL(8),DIMENSION(NI,NJ):: x, y
	   REAL(8) dx, dy, EPS, CFL, dt, A, Ur, Ul, Ut, Ub, Vt, Vb, Vr, Vl, Pl, Pr, Pt, Rt, Pb, MU, U0, Ru, Rv, Rp, MASS, IMPLX, IMPLY, Ue, Ve, Const, P0, R0, h_v, U1, U2
	   
	   INTEGER IU, IV, IP
	   
	   IU = 1
	   IV = 2
	   IP = 3
	   
	   open(IU, file = "Ru.txt")
	   open(IV, file = "Rv.txt")
	   open(IP, file = "Rp.txt")
	   
	   dt = CFL * (dx + dy) / (2.0 * U0)
	   
	   A = 1 / (U0**2)
	   
	   
	   Const = P0 / R0**(1.4)

		!Resp = 0
		!Resx = 0
		!Resy = 0
	   
	   do k = 1, NITER
		!Uold = U_c
		!Vold = V_c
		!Pold = P_c
		!Rold = R_c
		Ru = 0
		Rv = 0
		Rp = 0
		Resp = 0
		Resx = 0
		Resy = 0
		
		
		do j = 1, NJ - 1
			do i = 0, NI - 1 !?
				Ue = (U_c(i + 1, j) * R_c(i + 1, j) + U_c(i, j) * R_c(i, j)) * 0.5
				if (Ue .GE. 0) then
					Ut = U_c(i, j)
					Pt = P_c(i + 1, j)
					Vt = V_c(i, j)
					Rt = R_c(i, j)
				else
					Ut = U_c(i + 1, j)
					Pt = P_c(i, j)
					Vt = V_c(i + 1, j)
					Rt = R_c(i + 1, j)
				endif
				MASS = Ut * Rt
				IMPLX = Ue * Ut + Pt - 4.0 * MU * (U_c(i + 1, j) - U_c(i, j)) / (3.0 * dx)
				IMPLY = Ue * Vt - MU * (V_c(i + 1, j) - V_c(i, j)) / dx
				
				Resp(i, j) = Resp(i, j) + MASS / dx
				Resp(i + 1, j) = Resp(i + 1, j) - MASS / dx
				
				Resx(i, j) = Resx(i, j) + IMPLX / dx
				Resx(i + 1, j) = Resx(i + 1, j) - IMPLX / dx
				
				Resy(i, j) = Resy(i, j) + IMPLY / dx
				Resy(i + 1, j) = Resy(i + 1, j) - IMPLY / dx
				
				if (i > 0) then 
					Resx(i, j) = Resx(i, j) + 2.0 * MU * (V_c(i+1,j+1) - V_c(i-1,j+1) - V_c(i+1,j-1) + V_c(i-1,j-1)) / (12.0 * dx * dy)
					Resy(i, j) = Resy(i, j) - MU * (U_c(i+1,j+1) - U_c(i-1,j+1) - U_c(i+1,j-1) + U_c(i-1,j-1)) / (4.0 * dx * dy)
				endif
				
			enddo
		enddo
		
		do i = 1, NI - 1
			do j = 0, NJ - 1 !?
				Ve = (V_c(i, j + 1) * R_c(i, j + 1) + V_c(i, j) * R_c(i, j)) * 0.5
				if (Ve .GE. 0) then
					Ut = U_c(i, j)
					Pt = P_c(i, j + 1)
					Vt = V_c(i, j)
					Rt = R_c(i, j)
				else
					Ut = U_c(i, j + 1)
					Pt = P_c(i, j)
					Vt = V_c(i, j + 1)
					Rt = R_c(i, j + 1)
				endif
				MASS = Vt * Rt
				IMPLX = Ve * Ut - MU * (U_c(i, j + 1) - U_c(i, j)) / dy
				IMPLY = Ve * Vt + Pt - 4.0 * MU * (V_c(i, j + 1) - V_c(i, j)) / (3.0 * dy)
				
				Resp(i, j) = Resp(i, j) + MASS / dy
				Resp(i, j + 1) = Resp(i, j + 1) - MASS / dy
				
				Resx(i, j) = Resx(i, j) + IMPLX / dy
				Resx(i, j + 1) = Resx(i, j + 1) - IMPLX / dy
				
				Resy(i, j) = Resy(i, j) + IMPLY / dy
				Resy(i, j + 1) = Resy(i, j + 1) - IMPLY / dy
				
				if (j > 0) then 
					Resx(i, j) = Resx(i, j) - MU * (V_c(i+1,j+1) - V_c(i-1,j+1) - V_c(i+1,j-1) + V_c(i-1,j-1)) / (4.0 * dx * dy)
					Resy(i, j) = Resy(i, j) + 2.0 * MU * (U_c(i+1,j+1) - U_c(i-1,j+1) - U_c(i+1,j-1) + U_c(i-1,j-1)) / (12.0 * dx * dy)
				endif
			enddo
		enddo
		
		Resp = -dt * Resp * U0**2
		Resx = -dt * Resx
		Resy = -dt * Resy
		
	    P_c = P_c + Resp
		R_c = (P_c / Const)**(1.0/1.4)
		do i = 0, NI
			do j = 0, NJ
				U_c(i, j) = U_c(i, j) + Resx(i, j) / R_c(i, j)
				V_c(i, j) = V_c(i, j) + Resy(i, j) / R_c(i, j)
			enddo
		enddo
		
		
		
		!ГУ 
		
		!Нижняя граница:
		U_c(1:NI-1, 0) = U_c(1:NI-1, 1)
		V_c(1:NI-1, 0) = -V_c(1:NI-1, 1)
		P_c(1:NI-1, 0) = P_c(1:NI-1, 1)
		R_c(1:NI-1, 0) = R_c(1:NI-1, 1)
		
		!Входная граница:
		do j = 0, NJ
			if (Y_c(0, j) < h_V) then
				U_c(0, j) = U1
			else
				U_c(0, j) = U2
			endif
			V_c(0, j) = 0.0
			P_c(0, j) = P_c(1, j)
			R_c(0, j) = R_c(1, j)
		enddo
		
		!Верхнняя граница:
		U_c(1:NI-1, NJ) = U_c(1:NI-1, NJ - 1)
		V_c(1:NI-1, NJ) = -V_c(1:NI-1, NJ - 1)
		P_c(1:NI-1, NJ) = P_c(1:NI-1, NJ - 1)
		R_c(1:NI-1, NJ) = R_c(1:NI-1, NJ - 1)
		
		!Выходная граница:
		U_c(NI, :) = U_c(NI - 1, :)
		V_c(NI, :) = V_c(NI - 1, :)
		P_c(NI, :) = P0
		R_c(NI, :) = R0
		
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
		print*, k, Ru, Rv, Rp
	   enddo
	   
	   close(IU)
	   close(IV)
	   close(IP)

       END  SUBROUTINE

!************************************************************************************************
                    