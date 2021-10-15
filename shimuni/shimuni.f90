PROGRAM SHIMUNI
	IMPLICIT NONE
	CHARACTER(*), PARAMETER 				::	INPUT_FILE = "input.txt", &
												OUTPUT_FILE = "output.plt", &
												RESIDUALS_FILE = "residuals.dat"
	INTEGER(4), PARAMETER					::	COMMON_IO = 101
	REAL(8), DIMENSION(:), ALLOCATABLE 		::	NODE_ARR_X, &
												NODE_ARR_Y, &
												P_ARR
	REAL(8), DIMENSION(:, :), ALLOCATABLE	::	U_MATR, &
												V_MATR
	REAL(8)									::	START_X, &
												END_X, &
												START_Y, &
												END_Y, &
												DX, &
												DY, &
												MU, &
												U0, &
												V0, &
												P0, &
												EPS, &
												RE_X, &
												CF_ANALYT, &
												CF_THEOR
	INTEGER(4)								::	NUMB_X, &
												NUMB_Y, &
												NUMB_I, &
												NUMB_J, &
												I,		&
												J,		&
												SMAX
	
		OPEN(COMMON_IO, FILE=INPUT_FILE)
		CALL READ_INPUT_FILE(COMMON_IO, START_X, END_X, NUMB_X, &
						START_Y, END_Y, NUMB_Y, U0, V0, &
						P0, MU, EPS, SMAX)
		CLOSE(COMMON_IO)
		NUMB_I = NUMB_Y + 1
		NUMB_J = NUMB_X + 1
		
		ALLOCATE(NODE_ARR_X(NUMB_X))
		ALLOCATE(NODE_ARR_Y(NUMB_Y))
		ALLOCATE(U_MATR(NUMB_I, NUMB_J))
		ALLOCATE(V_MATR(NUMB_I, NUMB_J))
		ALLOCATE(P_ARR(NUMB_J))
		DX = CREATE_DIM(NODE_ARR_X, START_X, END_X, NUMB_X)
		DY = CREATE_DIM(NODE_ARR_Y, START_Y, END_Y, NUMB_Y)
		
		P_ARR(1) = P0

		DO I = 2, NUMB_Y, 1
			U_MATR(I, 1) = U0
			V_MATR(I, 1) = V0
		ENDDO
		U_MATR(1, 1) = -U_MATR(2, 1)
		V_MATR(1, 1) = -V_MATR(2, 1)
		U_MATR(NUMB_I, 1) = U_MATR(NUMB_Y, 1)
		V_MATR(NUMB_I, 1) = V_MATR(NUMB_Y, 1)
		
		OPEN(COMMON_IO, FILE=RESIDUALS_FILE)
		CALL SHIMUNI_SOLVER(COMMON_IO, NUMB_I, NUMB_J, SMAX, &
							DX, DY, MU, U0, V0, P0, &
							EPS, P_ARR, U_MATR, V_MATR)
		CLOSE(COMMON_IO)
		
		OPEN(COMMON_IO, FILE=OUTPUT_FILE)
		CALL WRITE_TO_TECPLOT(COMMON_IO, NUMB_X, NUMB_Y, &
						NODE_ARR_X, NODE_ARR_Y, P_ARR, &
						U_MATR, V_MATR)
		CLOSE(COMMON_IO)
		
		OPEN(COMMON_IO, FILE='cf.dat')
		DO J=2, NUMB_X, 1
			RE_X = U0 * (J - 0.5D0) * DX / MU
			CF_THEOR = 0.664D0 / SQRT(RE_X)
			CF_ANALYT = 4.0D0 * MU * U_MATR(2, J) / (DY * U0**2)
			WRITE(COMMON_IO, *) RE_X , CF_THEOR, CF_ANALYT
		ENDDO
		CLOSE(COMMON_IO)
		
		DEALLOCATE(P_ARR)
		DEALLOCATE(V_MATR)
		DEALLOCATE(U_MATR)
		DEALLOCATE(NODE_ARR_Y)
		DEALLOCATE(NODE_ARR_X)
		
		
	CONTAINS
		
		SUBROUTINE READ_INPUT_FILE(	IO, &
									START_X, &
									END_X, &
									NUMB_X, &
									START_Y, &
									END_Y, &
									NUMB_Y, &
									U0, &
									V0, &
									P0, &
									MU, &
									EPS, &
									SMAX)
									
		!READ THE INPUT FILE WITH INITIALIZATION PARAMETERS
		!IO - THE BUFFER ID
		!START_X, END_X, NUMB_X - A DEFINITION OF X-DIMENSION
		!START_Y, END_Y, NUMB_Y - A DEFINITION OF Y-DIMENSION
		!U0, V0, P0	- INPUT PARAMETERS
		!DENSITY - CONSTANT DENSITY OF FLUID
		!MU - CONSTANT DYNAMIC VISCOSITY COEFFICIENT
		!EPS - AN ACCEPTABLE MARGIN OF ERROR
		!SMAX - MAX NUMBER OF ITERATION PER VERTICAL LINE
		
			IMPLICIT NONE
			INTEGER(4), INTENT(IN)	::	IO
			INTEGER(4), INTENT(OUT)	::	NUMB_X, &
										NUMB_Y, &
										SMAX
			REAL(8), INTENT(OUT)	::	START_X, &
										END_X, &
										START_Y, &
										END_Y, &
										U0, &
										V0, &
										P0, &
										MU, &
										EPS
			READ(IO, *)	START_X, END_X, NUMB_X
			READ(IO, *) START_Y, END_Y, NUMB_Y
			READ(IO, *) U0, V0
			READ(IO, *) P0
			READ(IO, *) MU
			READ(IO, *) EPS
			READ(IO, *) SMAX
		
		END SUBROUTINE READ_INPUT_FILE
		
		
		
		SUBROUTINE TDMA(	A_ARR, &
							B_ARR, &
							C_ARR, &
							D_ARR, &
							X_ARR, &
							N)
							
		!SOLVES SYSTEM_MATRIX * X_ARR = D_ARR WHERE SYSTEM_MATRIX IS A 
		!TRIDIAGONAL MATRIX CONSISTING OF VERCTORS A_ARR, B_ARR, C_ARR
		!A_ARR - SUBDIAGONAL (MEANS IT IS THE DIAGONAL BELOW THE MAIN
		!		 DIAGONAL)
		!B_ARR - THE MAIN DIAGONAL
		!C_ARR - SUPERDIAGONAL (MEANS IT IS THE DIAGONAL ABOVE THE MAIN
		!		 DIAGONAL)
		!D_ARR - RIGHT PART OF THE SYSTEM
		!N - NUMBER OF EQUATIONS (LENGTH OF VECTORS)
		
			IMPLICIT NONE
			REAL(8), DIMENSION(:), INTENT(IN)	::	A_ARR, &
													B_ARR, &
													C_ARR, &
													D_ARR
			REAL(8), DIMENSION(:), INTENT(OUT)	::	X_ARR
			INTEGER(4), INTENT(IN)				::	N
			REAL(8), DIMENSION(N)				::	ALPHA_ARR, &
													BETA_ARR
			INTEGER(4)							::	I
			REAL(8)								::	DEN
						
			ALPHA_ARR(1) = -C_ARR(1) / B_ARR(1)
			BETA_ARR(1) = D_ARR(1) / B_ARR(1)
			DO I = 2, N - 1, 1
				DEN = B_ARR(I) + A_ARR(I) * ALPHA_ARR(I - 1)
				ALPHA_ARR(I) = -C_ARR(I) / DEN
				BETA_ARR(I) = (D_ARR(I) - A_ARR(I) * BETA_ARR(I - 1)) / DEN
			ENDDO
			
			DEN = B_ARR(N) + A_ARR(N) * ALPHA_ARR(N - 1)
			X_ARR(N) = (D_ARR(N) - A_ARR(N) * BETA_ARR(N - 1))/DEN
			
			DO I = N - 1, 1, -1
				X_ARR(I) = ALPHA_ARR(I) * X_ARR(I + 1) + BETA_ARR(I)
			ENDDO
			
		END SUBROUTINE TDMA
		
		
		
		REAL(8) FUNCTION TRAPEZOIDAL_INTEGRAL(	VAL_ARR, &
												NUMB, &
												DX)
		! NUMERICAL CALCULATE THE INTEGRAL ACCORDING TO THE TRAPEZOIDAL RULE
		! VAL_ARR - VALUES OF THE FUNCTION
		! NUMB - THE NUMBER OF CHUNKS
		! DX - COORDINATE STEP
			IMPLICIT NONE
			REAL(8), DIMENSION(:), INTENT(IN)	::	VAL_ARR
			INTEGER(4), INTENT(IN)				::	NUMB
			REAL(8), INTENT(IN)					::	DX
			INTEGER(4)							::	I
			
			TRAPEZOIDAL_INTEGRAL = (VAL_ARR(1)+ VAL_ARR(NUMB)) * 0.50D0
			DO I = 2, NUMB - 1, 1
				TRAPEZOIDAL_INTEGRAL = TRAPEZOIDAL_INTEGRAL + VAL_ARR(I)
			ENDDO
			TRAPEZOIDAL_INTEGRAL = TRAPEZOIDAL_INTEGRAL * DX
			
		END FUNCTION TRAPEZOIDAL_INTEGRAL
		


		REAL(8)  FUNCTION CREATE_DIM(	NODE_ARR_XI, &
										START_XI, &
										END_XI, &
										NUMB_XI)
		! NODE_ARR_XI - CREATE XI-COORDINATES
		! START_XI, END_XI - SECTION DEFINITION
		! NUMB_XI - NUMBER OF NODES
			IMPLICIT NONE
			REAL(8), DIMENSION(:), INTENT(OUT) 	:: 	NODE_ARR_XI
			INTEGER(4), INTENT(IN) 				:: 	NUMB_XI
			REAL(8), INTENT(IN) 				:: 	START_XI, &
													END_XI
			REAL(8)								::	DXI
			INTEGER(4)							::	I
			
			DXI = (END_XI - START_XI) / (NUMB_XI - 1)
			DO I = 1, NUMB_XI, 1
				NODE_ARR_XI(I) = START_XI + (I - 1) * DXI
			ENDDO
			CREATE_DIM = DXI
			
		END FUNCTION CREATE_DIM
		
		SUBROUTINE SHIMUNI_SOLVER(	IO, &
									NUMB_I, &
									NUMB_J, &
									SMAX, &
									DX, &
									DY, &
									MU, &
									U0, &
									V0, &
									P0, &
									EPS, &
									P_ARR, &
									U_MATR, &
									V_MATR)
		! SOLVE A SYSTEM OF PRANDTLS EQUATIONS
		! IO - THE BUFFER ID
		! NUMB_I, NUMB_J - DIMENSION OF SYSTEM
		! SMAX - ITERATION LIMITER
		! DX, DY - COORDIANTES STEPS
		! MU, - DYNAMIC VISCOSITY COEFFICIENT
		! U0, V0 - INPUT COMPONENTS OF VELOCITY
		! P0 - INPUT PRESSURE
		! EPS - AN ACCEPTABLE MARGIN OF ERROR
		! P_ARR - PRESSURE
		! U_MATR - X-VELOCITY
		! V_MATR - Y_VELOCITY
		
			IMPLICIT NONE
			INTEGER(4), INTENT(IN)						::	IO, &
															NUMB_I, &
															NUMB_J, &
															SMAX
			REAL(8), DIMENSION(:), INTENT(IN OUT)		::	P_ARR
			REAL(8), DIMENSION(:, :), INTENT(IN OUT)	::	U_MATR, &
															V_MATR
			REAL(8), INTENT(IN)							::	DX, &
															DY, &
															MU, &
															U0, &
															V0, &
															P0, &
															EPS
			REAL(8), DIMENSION(:), ALLOCATABLE			:: 	A_ARR, &
															B_ARR, &
															C_ARR, &
															D_ARR, &
															U_ARR, &
															V_ARR, &
															F_ARR, &
															Z_ARR, &
															W_ARR
			REAL(8)										::	G_NORM, &
															G_CURR, &
															G_IN, &
															MU_DIV, &
															W_INTEGRAL, &
															Z_INTEGRAL, &
															RES_U, &
															RES_V, &
															TEMP_RES, &
															V_NORM
			INTEGER(4)									::	I, &
															S, &
															K
															
			ALLOCATE(U_ARR(NUMB_I))
			ALLOCATE(V_ARR(NUMB_I))
			ALLOCATE(A_ARR(NUMB_I))
			ALLOCATE(B_ARR(NUMB_I))
			ALLOCATE(C_ARR(NUMB_I))
			ALLOCATE(D_ARR(NUMB_I))
			ALLOCATE(F_ARR(NUMB_I))
			ALLOCATE(Z_ARR(NUMB_I))
			ALLOCATE(W_ARR(NUMB_I))
			
			MU_DIV = MU / (DY * DY)
			
			DO I = 2, NUMB_I - 1, 1
				F_ARR(I) = -1.0D0 / DX
				U_ARR(I) = U_MATR(I, 1)
				V_ARR(I) = V_MATR(I, 1)
			ENDDO
			F_ARR(1) = 0.0D0
			F_ARR(NUMB_I) = 0.0D0
			U_ARR(1) = -U_ARR(2)
			U_ARR(NUMB_I) = U_ARR(NUMB_I - 1)
			V_ARR(1) = -V_ARR(2)
			V_ARR(NUMB_I) = V_ARR(NUMB_I - 1)
			V_NORM = SQRT(U0 * U0 + V0 * V0)
			
			G_IN = TRAPEZOIDAL_INTEGRAL(U_ARR, NUMB_I, DY)
			G_NORM = G_IN
			K = 0
			DO J = 2, NUMB_J, 1
				P_ARR(J) = P_ARR(J - 1)
				S = 0
				RES_U = 2.0D0 * EPS
				RES_V = 2.0D0 * EPS
				DO WHILE ((S < SMAX).AND.((RES_U.GE.EPS).OR.(RES_V.GE.EPS)))
					K = K + 1
					S = S + 1
					A_ARR(1) = 0.0D0
					B_ARR(1) = 1.0D0
					C_ARR(1) = 1.0D0
					D_ARR(1) = 0.0D0
					
					DO I = 2, NUMB_I-1, 1
						A_ARR(I) = -0.5D0 * V_ARR(I - 1) / DY - MU_DIV
						B_ARR(I) = U_ARR(I) / DX + 2.0D0 * MU_DIV
						C_ARR(I) = 0.5D0 * V_ARR(I + 1) / DY - MU_DIV
						D_ARR(I) = (U_MATR(I, J - 1) ** 2.0D0 + P_ARR(J - 1)) / DX
					ENDDO
					
					A_ARR(NUMB_I) = -1.0D0
					B_ARR(NUMB_I) = 1.0D0
					C_ARR(NUMB_I) = 0.0D0
					D_ARR(NUMB_I) = 0.0D0
					
					!OUTER FLOW CONDITION
					!A_ARR(NUMB_I) = 0.0D0
					!B_ARR(NUMB_I) = 1.0D0
					!C_ARR(NUMB_I) = 0.0D0
					!D_ARR(NUMB_I) = U0
					
					CALL TDMA(A_ARR, B_ARR, C_ARR, D_ARR, Z_ARR, NUMB_I)
					CALL TDMA(A_ARR, B_ARR, C_ARR, F_ARR, W_ARR, NUMB_I)
					W_INTEGRAL = TRAPEZOIDAL_INTEGRAL(W_ARR, NUMB_I, DY)
					Z_INTEGRAL = TRAPEZOIDAL_INTEGRAL(Z_ARR, NUMB_I, DY)
					P_ARR(J) = (G_IN + (V_ARR(1) - V_ARR(NUMB_I)) * DX - Z_INTEGRAL) / W_INTEGRAL
					RES_U = 0.0D0
					DO I = 1, NUMB_I, 1
						TEMP_RES = U_ARR(I)
						U_ARR(I) = Z_ARR(I) + W_ARR(I) * P_ARR(J)
						RES_U = MAX(RES_U, ABS(TEMP_RES - U_ARR(I))/V_NORM)
					ENDDO
					
					G_CURR = TRAPEZOIDAL_INTEGRAL(U_ARR, NUMB_I, DY)
					V_ARR(1) = 0.25D0*DY *(U_ARR(2) + U_ARR(1) - U_MATR(2,J-1) - U_MATR(1,J-1))/DX
					RES_V = V_ARR(1)
					DO I = 2, NUMB_I, 1
						TEMP_RES = V_ARR(I)
						V_ARR(I) = V_ARR(I - 1) - 0.5D0*DY * (U_ARR(I) + U_ARR(I-1) - U_MATR(I,J-1) - U_MATR(I-1,J-1))/DX
						RES_V = MAX(RES_V, ABS(TEMP_RES - V_ARR(I))/V_NORM)
					ENDDO
					
					TEMP_RES = G_CURR/G_NORM
					WRITE(*, *) K, S, (J - 0.5D0)* DX, TEMP_RES, RES_U, RES_V
					WRITE(IO, *) K, S, (J - 0.5D0)* DX, TEMP_RES, RES_U, RES_V
				ENDDO
				
				G_IN = G_CURR
				DO I = 1, NUMB_I, 1
					U_MATR(I, J) = U_ARR(I)
					V_MATR(I, J) = V_ARR(I)
				ENDDO	
			ENDDO
			
			DEALLOCATE(W_ARR)
			DEALLOCATE(Z_ARR)
			DEALLOCATE(F_ARR)
			DEALLOCATE(D_ARR)
			DEALLOCATE(C_ARR)
			DEALLOCATE(B_ARR)
			DEALLOCATE(A_ARR)
			DEALLOCATE(V_ARR)
			DEALLOCATE(U_ARR)
															
		END SUBROUTINE SHIMUNI_SOLVER
		
		SUBROUTINE WRITE_TO_TECPLOT(IO, 		&
									NUMB_X,		&
									NUMB_Y,		&
									NODE_ARR_X,	&
									NODE_ARR_Y,	&
									P_ARR,		&
									U_MATR,		&
									V_MATR)
		! OUTPUT DATA TO TECPLOT IN A FINITE VOLUME FORMAT
		! IO - THE BUFFER ID
		! NUMB_X, NUMB_Y - NUMBER OF NODES (NUMB_Y x NUMB_X)
		! NODE_ARR_X - X-COORDINATES OF NODES
		! NODE_ARR_Y - Y-COORDINATES OF NODES
		! P_ARR - PRESSURE
		! U_MATR - X-COMPONENT OF THE VELOCITY VECTOR
		! V_MATR - Y-COMPONENT OF THE VELOCITY VECTOR
			INTEGER(4), INTENT(IN)				::	IO
			INTEGER(4), INTENT(IN)				::	NUMB_X, &
													NUMB_Y
			REAL(8), DIMENSION(:), INTENT(IN)	::	NODE_ARR_X, &
													NODE_ARR_Y, &
													P_ARR
			REAL(8), DIMENSION(:, :), INTENT(IN)::	U_MATR, &
													V_MATR
			INTEGER(4)							::	I, J
													
			WRITE(IO, *)	'VARIABLES = "X", "Y", "P", "U", "V"'
			WRITE(IO, *)	'ZONE I=', NUMB_Y, ', J=', NUMB_X, ', DATAPACKING=BLOCK, VARLOCATION=([3-5]=CELLCENTERED)'
			WRITE(IO, '(100F14.7)')	((NODE_ARR_X(J), I=1, NUMB_Y), J=1, NUMB_X)
			WRITE(IO, '(100F14.7)')	((NODE_ARR_Y(I), I=1, NUMB_Y), J=1, NUMB_X)
			WRITE(IO, '(100F14.7)') ((P_ARR(J), I=2, NUMB_Y), J=2, NUMB_X)
			WRITE(IO, '(100F14.7)') U_MATR(2:NUMB_Y, 2:NUMB_X)
			WRITE(IO, '(100F14.7)') V_MATR(2:NUMB_Y, 2:NUMB_X)
		END SUBROUTINE WRITE_TO_TECPLOT
		
END PROGRAM SHIMUNI