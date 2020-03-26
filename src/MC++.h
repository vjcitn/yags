/* header file for MC++ library */

/* /usr13/stdevs/stdev0f/MC+.95.2/SCCS/s.MC++.h */
/* MC++.h (c) V. Carey, ver. @(#) 1.7 95/10/16 */

#ifndef HAS_MCPP_H
#define HAS_MCPP_H

#include <stdlib.h>
#include <setjmp.h>
#include <math.h>

#define FOR_S
#define assigned_el( matrix , i , j ) *( matrix . mathead() + ( i * matrix . cols() ) + j )
#define put_el( matrix , i , j ) *( matrix . mathead() + ( i * matrix . cols() ) + j )
#define set_el( matrix , i , j ) *( matrix . mathead() + ( i * matrix . cols() ) + j )

#ifndef FOR_S
#include <iostream.h>

	/* this macro sets up an output stream -- */
	/* arg1 : char* (filename) */
	/* arg2 : stream name (undeclared token) */
	/* arg3 : random name for buffer -- needs to be unique */
	/* after running it, you can use streamname << ... */
#define set_ostream( filename , stream , BUF ) \
	filebuf BUF ; \
	if  ( BUF.open( filename , output ) == 0 ) \
		{ \
		cerr << "MC++: set_ostream: can't open " << filename << "\nDies.\n"; \
		exit(1); \
		} \
\
	ostream stream( &BUF );

#else
#include <stdio.h>
#endif


#ifndef FOR_S
#define PRINT cout.flush() ; cout <<   /* helps a little with synchrony */
#else
#define PRINT(x)
#endif

extern jmp_buf mcpp_env;
#define error_signal( buf, code ) longjmp( buf, code )

#define  MSORT_FAIL_NO_COL 1
#define  CHOL_FAIL_NOT_SQR 2
#define CHOL_FAIL_ILLG_DIAG 3
#define CHOL_FAIL_NON_POS 4
#define  CLUSCNT_FAIL_SPLIT_ERR 5
#define  COMMON_FAIL_DIM_ERR 6
#define  COMMON_FAIL_NEED_VECT 7
#define  CONCAT_FAIL_ROW_AGRMNT 8
#define  CONCAT_FAIL_COL_AGRMNT 9
#define  COVLAG_FAIL_MATX_SIZE 10
#define  DELCOL_FAIL_NO_COL 11
#define  DELROW_FAIL_NO_ROW 12
#define  ELWISE_FAIL_DIM_AGRMNT 13
#define  EXTRACT_FAIL_NOT_SQR 14
#define  FORM_FAIL_NOT_VECT 15
#define  FSWEEP_FAIL_TOO_SMALL 16
#define  FSWEEP_FAIL_NOT_SEMIDEF 17
#define  FSWEEPINP_FAIL_TOO_SMALL 18
#define  FSWEEPINP_FAIL_NOT_SEMIDEF 19
#define  FSWEEPINP_FAIL_TRIANGLE 20
#define  LULINP_FAIL_NOT_SQR 21
#define  MATTIMES_FAIL_NOT_COLVEC 22
#define  MATTIMES_FAIL_DIM_AGRMNT 23
#define  MATMATH_FAIL_NEG_MATRIX 24
#define  MATPOW_FAIL_NOT_SQR 25
#define  MATPOW_FAIL_INVERSE_ONLY 26
#define  MAXREL_FAIL_NONMATCH_ARGS 27
#define  MINUS_FAIL_DIM_AGRMNT 28
#define  MINUSEQ_FAIL_DIM_AGRMNT 29
#define  MULT_FAIL_DIM_AGRMNT 30
#define  PLUG_FAIL_DIM_AGRMNT 31
#define  PLUS_FAIL_DIM_AGRMNT 32
#define  PLUSEQ_FAIL_DIM_AGRMNT 33
#define  ROWMUL_FAIL_NOT_COLVEC 34
#define  ROWMUL_FAIL_DIM_AGRMT 35
#define  ROWSEG_FAIL_DIM_AGRMT 36
#define  SEQ_FAIL_INFINITE_SEQ 37
#define  SPLIT_FAIL_NOT_COLVEC 38
#define  SPLIT_FAIL_DIM_AGRMNT 39
#define  SPLIT_FAIL_NO_POINTERS 40
#define  SUBMAT_FAIL_NOT_ROWVEC 41
#define  SUBMAT_FAIL_NOT_COLVEC 42
#define  SUBMAT_FAIL_CHECK_ROWS 43
#define  SUBMAT_FAIL_CHECK_COLS 44
#define  TOEPLITZ_FAIL_DIM_AGRMNT 45
#define  VEC2DIAG_FAIL_NEED_DIM1 46
#define  MATREAD_FAIL_CANT_OPEN 47
#define  MATREAD_FAIL_NOT_RECT 48
#define  MATMATH_FAIL_NEG_SQRT 49
#define SUPPLIED_S_TOO_LARGE 101
#define TOO_MANY_DISTINCT_S 102
#define NON_CONVERGENCE 103
#define INVALID_LOGIT_ARG 104
#define INVALID_ANTILOGIT_ARG 105

#endif

/* end of header MC++.h */

