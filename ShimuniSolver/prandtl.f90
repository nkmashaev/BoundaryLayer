PROGRAM PRANDTL
	USE GRID_MODULE
	USE LINALG_MODULE
	USE TECPLOT_MODULE
	USE CALCULUS_MODULE
	IMPLICIT NONE
	
	INTEGER(4), PARAMETER					::	IO=101, IO_TEMP=102
	REAL(8), PARAMETER						:: 	EPS=1.E-7
	REAL(8), DIMENSION(:), ALLOCATABLE 		::	NODE_ARR_X, &
												NODE_ARR_Y, &
												A_ARR, &
												B_ARR, &
												C_ARR, &
												D_ARR, &
												P_ARR, &
												U_ARR, &
												V_ARR, &
												F_ARR, &
												Z_ARR, &
												W_ARR
	REAL(8), DIMENSION(:, :), ALLOCATABLE	::	U_MATR, &
												V_MATR
	REAL(8)									::	START_X, &
												END_X, &
												START_Y, &
												END_Y, &
												DX, &
												DY, &
												MU, &
												MU_DIV, &
												G_IN, &
												G_CURR, &
												G_NORM, &
												W_INTEGRAL, &
												Z_INTEGRAL, &
												RES, &
												MAX_U_RES, &
												MAX_V_RES, &
												U0, &
												RE_X, &
												CF_AN, &
												CF_NUM
	INTEGER(4)								::	NUMB_X, &
												NUMB_Y, &
												NUMB_I, &
												NUMB_J, &
												I,		&
												J,		&
												K,		&
												S,		&
												BACKUP,	&
												SMAX,	&
												JSTART, &
												KSTART
		MU = 1.0E-5
		U0 = 0.01D0
		BACKUP = 10000
		G_NORM = -1.0D0
		
		SMAX = 500
		START_X = 0.0D0
		END_X = 0.1D0
		START_Y = 0.0D0
		END_Y = 0.1D0
		NUMB_X = 201
		NUMB_Y = 201
		NUMB_I = NUMB_Y + 1
		NUMB_J = NUMB_X + 1
		JSTART = 2
		KSTART = 0
		ALLOCATE(NODE_ARR_X(NUMB_X))
		ALLOCATE(NODE_ARR_Y(NUMB_Y))
		ALLOCATE(U_MATR(NUMB_I, NUMB_J))
		ALLOCATE(V_MATR(NUMB_I, NUMB_J))
		ALLOCATE(P_ARR(NUMB_J))
		DX = CREATE_DIM(NODE_ARR_X, START_X, END_X, NUMB_X)
		DY = CREATE_DIM(NODE_ARR_Y, START_Y, END_Y, NUMB_Y)
		U_MATR(1:NUMB_I, 1:NUMB_J) = U0
		U_MATR(1, 1:NUMB_J) = -U_MATR(2, 1:NUMB_J)
		V_MATR(1:NUMB_I, 1:NUMB_J) = 0.0D0
		V_MATR(1, 1:NUMB_J) = -V_MATR(2, 1:NUMB_J)
		P_ARR(1:NUMB_J) = 0.0D0
		
		!OPEN(IO, FILE='temp_output.dat')
		!CALL READ_TEMP_SOLUTION(IO, &
		!						KSTART, &
		!						JSTART, &
		!						NUMB_X, &
		!						NUMB_Y, &
		!						NUMB_I, &
		!						NUMB_J, &
		!						DX, &
		!						DY, &
		!						G_NORM, &
		!						NODE_ARR_X, &
		!						NODE_ARR_Y, &
		!						P_ARR, &
		!						U_MATR, &
		!						V_MATR)
		!CLOSE(IO)
		
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
		F_ARR(1) = 0.0D0
		F_ARR(2:NUMB_Y) = -1.0D0 / DX
		F_ARR(NUMB_I) = 0.0D0
		U_ARR(1 : NUMB_I) = U_MATR(1 : NUMB_I, JSTART - 1)
		V_ARR(1 : NUMB_I) = V_MATR(1 : NUMB_I, JSTART - 1)
		
		
		G_IN = TRAPEZOIDAL_INTEGRAL(U_ARR, NUMB_I, DY)
		IF (G_NORM.LT.0.0D0) THEN
			G_NORM = G_IN
		ENDIF
		WRITE(*, *) G_IN
		
		W_INTEGRAL = 0.0D0
		Z_INTEGRAL = 0.0D0
		
		OPEN(IO, FILE='output.plt')
		CALL FV_OUTPUT(IO, NUMB_X, NUMB_Y, NODE_ARR_X, & 
						NODE_ARR_Y, P_ARR, U_MATR, V_MATR)
		CLOSE(IO)
		K = KSTART
		OPEN(IO, FILE='residuals.dat')
		DO J = JSTART, NUMB_J, 1
			P_ARR(J) = P_ARR(J - 1)
			S = 0
			DO WHILE(S < SMAX)
				K = K + 1
				S = S + 1
				!symmetry
				!A_ARR(1) = 0.0D0
				!B_ARR(1) = 1.0D0
				!C_ARR(1) = -1.0D0
				!D_ARR(1) = 0.0D0
			
				!wall
				A_ARR(1) = 0.0D0
				B_ARR(1) = 1.0D0
				C_ARR(1) = 1.0D0
				D_ARR(1) = 0.0D0
				
				DO I = 2, NUMB_Y, 1
					A_ARR(I) = -0.5D0 * V_ARR(I - 1) / DY - MU_DIV
					B_ARR(I) = U_ARR(I) / DX + 2.0D0 * MU_DIV
					C_ARR(I) = 0.5D0 * V_ARR(I + 1) / DY - MU_DIV
					D_ARR(I) = (U_MATR(I, J-1) * U_MATR(I, J-1) + P_ARR(J - 1)) / DX
				ENDDO
				
				!A_ARR(NUMB_I) = 0.0D0
				!B_ARR(NUMB_I) = U0
				!C_ARR(NUMB_I) = 0.0D0
				!D_ARR(NUMB_I) = 0.0D0				
				
				!symmetry
				A_ARR(NUMB_I) = -1.0D0
				B_ARR(NUMB_I) = 1.0D0
				C_ARR(NUMB_I) = 0.0D0
				D_ARR(NUMB_I) = 0.0D0
				
				!wall
				!A_ARR(NUMB_I) = 1.0D0
				!B_ARR(NUMB_I) = 1.0D0
				!C_ARR(NUMB_I) = 0.0D0
				!D_ARR(NUMB_I) = 0.0D0
			
				CALL TDMA(A_ARR, B_ARR, C_ARR, D_ARR, Z_ARR, NUMB_I)
				CALL TDMA(A_ARR, B_ARR, C_ARR, F_ARR, W_ARR, NUMB_I)
				W_INTEGRAL = TRAPEZOIDAL_INTEGRAL(W_ARR, NUMB_I, DY)
				Z_INTEGRAL = TRAPEZOIDAL_INTEGRAL(Z_ARR, NUMB_I, DY)
				P_ARR(J) = (G_IN + (V_ARR(1) - V_ARR(NUMB_I)) * DX - Z_INTEGRAL) / W_INTEGRAL
				MAX_U_RES = 0.0D0
				DO I = 1, NUMB_I, 1
					RES = U_ARR(I)
					U_ARR(I) = Z_ARR(I) + W_ARR(I) * P_ARR(J)
					MAX_U_RES = MAX(MAX_U_RES, ABS(RES - U_ARR(I))/MAX(ABS(RES), EPS))
				ENDDO
				G_CURR = TRAPEZOIDAL_INTEGRAL(U_ARR, NUMB_I, DY)
				V_ARR(1) = 0.25D0*DY *(U_ARR(2) + U_ARR(1) - U_MATR(2,J-1) - U_MATR(1,J-1))/DX
				MAX_V_RES = 0.0D0
				DO I = 2, NUMB_I, 1
					RES = V_ARR(I)
					V_ARR(I) = V_ARR(I - 1) - 0.5D0*DY * (U_ARR(I) + U_ARR(I-1) - U_MATR(I,J-1) - U_MATR(I-1,J-1))/DX
					MAX_V_RES = MAX(MAX_V_RES, ABS(RES - V_ARR(I))/MAX(ABS(RES), EPS))
				ENDDO

				WRITE(*, *) K, S, NODE_ARR_X(J-1) + 0.5D0* DX, G_CURR/G_NORM, MAX_U_RES, MAX_V_RES
				WRITE(IO, *) K, S, NODE_ARR_X(J-1) + 0.5D0* DX, G_CURR/G_NORM, MAX_U_RES, MAX_V_RES
				IF (MOD(K, BACKUP).EQ.0) THEN
					WRITE(*, *) 'DOING INTERMEDIATE DATA STORAGE...'
					OPEN(IO_TEMP, FILE='temp_output.dat')
					CALL SAVE_TEMP_SOLUTION(IO_TEMP, K-S, J, NUMB_X, NUMB_Y, &
											DX, DY, G_NORM, NODE_ARR_X, &
											NODE_ARR_Y, P_ARR, U_MATR, V_MATR)
					CLOSE(IO_TEMP)
					
					OPEN(IO_TEMP, FILE='output.plt')
					CALL FV_OUTPUT(IO_TEMP, NUMB_X, NUMB_Y, &
									NODE_ARR_X, NODE_ARR_Y, &
									P_ARR, U_MATR, V_MATR)
					CLOSE(IO_TEMP)
					
				ENDIF
			ENDDO
			G_IN = G_CURR
			DO I = 1, NUMB_I, 1
				U_MATR(I, J) = U_ARR(I)
				V_MATR(I, J) = V_ARR(I)
			ENDDO	
		ENDDO
		
		CLOSE(IO)
		
		WRITE(*, *) 'DOING FINAL DATA SAVING...'
		OPEN(IO, FILE='output.plt')
		CALL FV_OUTPUT(IO, NUMB_X, NUMB_Y, NODE_ARR_X, & 
						NODE_ARR_Y, P_ARR, U_MATR, V_MATR)
		CLOSE(IO)
		
		OPEN(IO, FILE='cf.dat') 
		DO J = 2, NUMB_X, 1
			RE_X = (NODE_ARR_X(J) - DX * 0.50D0) * U0 / MU
			CF_AN = 0.664D0 / SQRT(RE_X)
			CF_NUM = 2.0D0 * MU * (U_MATR(2, J) - U_MATR(1, J))/(DY * U0 * U0)
			WRITE(IO, '(4F14.7)') NODE_ARR_X(J) - DX * 0.50D0, RE_X, CF_AN, CF_NUM
		ENDDO
		CLOSE(IO)
		DEALLOCATE(W_ARR)
		DEALLOCATE(Z_ARR)
		DEALLOCATE(F_ARR)
		DEALLOCATE(D_ARR)
		DEALLOCATE(C_ARR)
		DEALLOCATE(B_ARR)
		DEALLOCATE(A_ARR)
		DEALLOCATE(V_ARR)
		DEALLOCATE(U_ARR)
		DEALLOCATE(P_ARR)
		DEALLOCATE(V_MATR)		
		DEALLOCATE(U_MATR)
		DEALLOCATE(NODE_ARR_Y)
		DEALLOCATE(NODE_ARR_X)
		
	CONTAINS
	
		SUBROUTINE SAVE_TEMP_SOLUTION(IO, &
									ITER, &
									CURR_J, &
									NUMB_X,		&
									NUMB_Y,		&
									DX,	&
									DY, &
									G_NORM, &
									NODE_ARR_X,	&
									NODE_ARR_Y,	&
									P_ARR,		&
									U_MATR,		&
									V_MATR)
			INTEGER(4), INTENT(IN)				::	IO, &
													ITER, &
													CURR_J, &
													NUMB_X, &
													NUMB_Y
			REAL(8), DIMENSION(:), INTENT(IN)	::	NODE_ARR_X, &
													NODE_ARR_Y, &
													P_ARR
			REAL(8), DIMENSION(:, :), INTENT(IN)::	U_MATR, &
													V_MATR
			REAL(8), INTENT(IN)					::	DX, &
													DY, &
													G_NORM
			INTEGER(4)							::	I, &
													J
			WRITE(IO, *) ITER
			WRITE(IO, *) CURR_J
			WRITE(IO, *) NUMB_X, NUMB_Y
			WRITE(IO, *) DX, DY
			WRITE(IO, *) G_NORM
			WRITE(IO, '(F14.7)') (NODE_ARR_X(J), J=1, NUMB_X)
			WRITE(IO, '(F14.7)') (NODE_ARR_Y(I), I=1, NUMB_Y)
			WRITE(IO, '(F14.7)') (P_ARR(J), J=1, NUMB_X + 1)
			WRITE(IO, '(F14.7)') U_MATR(1:NUMB_Y+1, 1:NUMB_X+1)
			WRITE(IO, '(F14.7)') V_MATR(1:NUMB_Y+1, 1:NUMB_X+1)
		END SUBROUTINE SAVE_TEMP_SOLUTION
		
		SUBROUTINE READ_TEMP_SOLUTION(IO, &
									KSTART, &
									JSTART, &
									NUMB_X, &
									NUMB_Y, &
									NUMB_I, &
									NUMB_J, &
									DX, &
									DY, &
									G_NORM, &
									NODE_ARR_X, &
									NODE_ARR_Y, &
									P_ARR, &
									U_MATR, &
									V_MATR)
			INTEGER(4), INTENT(IN)					::	IO
			INTEGER(4), INTENT(OUT)					::	KSTART, &
														JSTART, &
														NUMB_X, &
														NUMB_Y, &
														NUMB_I, &
														NUMB_J
			REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(OUT)	::	NODE_ARR_X, &
														NODE_ARR_Y, &
														P_ARR
			REAL(8), DIMENSION(:, :), ALLOCATABLE, INTENT(OUT)	::	U_MATR, &
														V_MATR
			REAL(8), INTENT(OUT)					::	DX,	& 
														DY, &
														G_NORM
			INTEGER(4)								::	I,	& 
														J
			
			READ(IO, *)	KSTART
			READ(IO, *) JSTART
			READ(IO, *) NUMB_X, NUMB_Y
			NUMB_I = NUMB_Y + 1
			NUMB_J = NUMB_X + 1
			READ(IO, *) DX, DY
			READ(IO, *) G_NORM
			ALLOCATE(NODE_ARR_X(NUMB_X))
			ALLOCATE(NODE_ARR_Y(NUMB_Y))
			ALLOCATE(U_MATR(NUMB_I, NUMB_J))
			ALLOCATE(V_MATR(NUMB_I, NUMB_J))
			ALLOCATE(P_ARR(NUMB_J))
			
			READ(IO, *) (NODE_ARR_X(J), J=1, NUMB_X)
			READ(IO, *) (NODE_ARR_Y(I), I=1, NUMB_Y)
			READ(IO, *) (P_ARR(J), J=1, NUMB_J)
			READ(IO, *) ((U_MATR(I, J), I=1, NUMB_I), J=1, NUMB_J)
			READ(IO, *) ((V_MATR(I, J), I=1, NUMB_I), J=1, NUMB_J)
		END SUBROUTINE READ_TEMP_SOLUTION

END PROGRAM PRANDTL