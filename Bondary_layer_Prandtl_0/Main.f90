    Program Pr
        Implicit none
		
 
        INTEGER, PARAMETER:: IO = 12 ! input-output unit
        INTEGER NI, NJ, NITER, SMAX
        INTEGER I,J
        REAL L,H,U0,MU,Nu,R0,P0
        REAL dx,dy,CFL,EPS
        REAL,ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
        REAL,ALLOCATABLE :: U_n(:,:),V_n(:,:),P_n(:,:),R_n(:,:)
		REAL, ALLOCATABLE :: Cf(:), Cf_B(:)

        write(*,*) 'Read input file'
        open(IO,FILE='Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) EPS
        read(IO,*) SMAX

        read(IO,*) U0
        read(IO,*) MU
        read(IO,*) R0
        read(IO,*) P0
        CLOSE(IO)
   
        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates

!----------------- Node variables -----------------------------
        allocate(U_n(NI,NJ))  ! Velocity U
        allocate(V_n(NI,NJ))  ! Velocity V
        allocate(P_n(NI,NJ))  ! Pressure
		
		allocate(Cf(NI), Cf_B(NI))

!----------------- Coordinate of nodes ------------------------
        dx=L/(NI-1)
        dy=H/(NJ-1)

        DO I=1,NI
          DO J=1,NJ
            X_Node(I,J)=(I-1)*dx
            Y_Node(I,J)=(J-1)*dy
          END DO
        END DO

!----------------- Parameters ------------------------

        NU=MU/R0

        write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
        write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
        write(*,*)'ReL= ', U0*L/NU
        pause    

!----------------- Initial fields -----------------------------

        DO I=1,NI
          DO J=1,NJ
            U_n(I,J)=U0
            V_n(I,J)=1.0e-5
            P_n(I,J)=P0
          ENDDO
        ENDDO
		

!---------------- Solve Prandtl equations ---------------------
 
        write(*,*) 'Solve Prandtl equations'      
        call  Prandtl(NI, NJ, X_Node, Y_Node, U_n, V_n, P_n, U0, EPS, SMAX, dx, dy, NU)
		call Find_tw(NI, NJ, X_Node, Y_Node, Cf, U_n, dy, MU)
		
		DO I=1,NI
          Cf_B(I) = 0.664 / sqrt(U0 * X_Node(I, 1) / NU)
        ENDDO
		
		Cf = 2 * Cf / (R0 * U0**2)

 !----------------- Output data ------------------------------
 
        write(*,*) 'Output data' 
        Open(IO,FILE='data_pr.plt')
        Call Output_Fields(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P_n)
		
		call Output_Cf(NI, NJ, X_Node * U0 / NU, Cf, Cf_B)
		
        Close(IO)
     
    END PROGRAM
	
	!include 'Prandtl.f90'
	!include 'Output.f90'
	!include 'Functions.f90'
   
