SUBMODULE (LINALG_MODULE) LINALG_LYNSYS

	CONTAINS
	
		MODULE SUBROUTINE TDMA( A_ARR, &
								B_ARR, &
								C_ARR, &
								D_ARR, &
								X_ARR, &
								N)
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

END SUBMODULE LINALG_LYNSYS