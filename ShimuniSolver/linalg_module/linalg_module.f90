MODULE LINALG_MODULE
	IMPLICIT NONE
	
	INTERFACE

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
			!SOLVES SYSTEM_MATRIX * X_ARR = D_ARR WHERE SYSTEM_MATRIX IS A 
			!TRIDIAGONAL MATRIX CONSISTING OF VERCTORS A_ARR, B_ARR, C_ARR
			!A_ARR - SUBDIAGONAL (MEANS IT IS THE DIAGONAL BELOW THE MAIN
			!		 DIAGONAL)
			!B_ARR - THE MAIN DIAGONAL
			!C_ARR - SUPERDIAGONAL (MEANS IT IS THE DIAGONAL ABOVE THE MAIN
			!		 DIAGONAL)
			!D_ARR - RIGHT PART OF THE SYSTEM
			!N - NUMBER OF EQUATIONS (LENGTH OF VECTORS)
		END SUBROUTINE TDMA

	END INTERFACE
	
END MODULE LINALG_MODULE