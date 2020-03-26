#include <stdio.h>
#include <string.h>
#include "MC++.h"

/* #define EXPERIMENTAL_USERCOR */

jmp_buf mcpp_env;
 
#define GetUjStruc( dBtag, Jtag, Btag, invfn, alpmat, niint, timmat, eps ) \
 matrix dBtag = diffL ( invfn , alpmat , niint , timmat, eps ); \
 matrix Jtag = form_diag ( extract_diag ( LDLt ( invfn ( alpmat, niint, timmat ) ) ) ); \
 matrix Btag = lowutri ( LDLt ( invfn ( alpmat, niint, timmat ) ) ) ;

#define GetUjKStruc( dBtag, Jtag, Btag, invfn, alpmat, niint, timmat, eps ) \
 matrix dBtag = diffLK ( invfn , alpmat , niint , timmat, eps ); \
 matrix Jtag = form_diag ( extract_diag ( LDLt (K* invfn ( alpmat, niint, timmat )*K ) ) ); \
 matrix Btag = lowutri ( LDLt ( K*invfn ( alpmat, niint, timmat )*K ) ) ;


 class matrix {
 
 		int nrows, ncols;
 		double *data;
 	public:
 		matrix();
 		matrix(const matrix&); // for initialztn of uninit vbl p. 180
 		~matrix();  
 
 		friend matrix newmat(int i, int j);
 
 		double el(int i, int j) { return *(data+(ncols*i)+j); };
 		matrix submat( matrix row_req , matrix col_req );
 		matrix submat( matrix row_req , int dummy  );
 		matrix submat( int dummy , matrix col_req ); 
 
 		double *mathead() { return data; };
 
 		int rows() { return nrows; };
 		int cols() { return ncols; };
 
 		friend matrix operator*( matrix , matrix );
 		friend matrix operator*( double , matrix );
 		friend matrix operator*( matrix , double );
 
 		friend matrix operator+( matrix , matrix );
 		friend void operator+=( matrix& , matrix& );
 		friend matrix operator-( matrix& , matrix& );
 		friend void operator-=( matrix& , matrix& );
 
 		friend matrix operator||( matrix , matrix );  /* H-CONCAT */
 		friend matrix operator/( matrix& , matrix& );   /* V-CONCAT */
 
 		friend matrix matlog( matrix& );
 		friend matrix matexp( matrix& );
 		friend matrix matsqrt( matrix& );
 		friend matrix matabs( matrix& );
 		friend matrix chol( matrix );
 		friend double elsum( matrix );
 
 		void operator=( matrix );  /* assignment */
 		friend void kill( matrix& );
 
 		};
 
 matrix::matrix()     /* basic constructor */
 	{            /* initialize metadata, prepare for assn or newmat */
 	nrows = 0;
 	ncols = 0;
 	};
 
 matrix::matrix(const matrix& inmat)   /* constructor for initialization by existing object */
 	{
 	nrows = inmat.nrows;
 	ncols = inmat.ncols;
 	data = new double[ nrows*ncols ];
 	for ( int i = 0 ; i < nrows*ncols ; i++ )
 		{
 		data[i] = inmat.data[i];
 		}
 	};
 
 matrix::~matrix()      /* destructor */
 	{
 	if ( nrows*ncols > 0 ) delete data;
 	};
 
 void matrix::operator=( matrix inmat )  // assignment !!
 	{				 // assumption -- this refers to lhs
 	//if ( *this == inmat ) return;    // inmat=inmat
 	if ( this->nrows > 0 ) delete data; /* cleans old mat */
 	this->nrows = inmat.nrows;
 	this->ncols = inmat.ncols;
 	int nelem = inmat.nrows*inmat.ncols ;
 	this->data = new double[ nelem ];
 	for ( int i = 0 ; i < nelem ; i++ )
 		{
 		this->data[i] = inmat.data[i];
 		}
 	};
 
 matrix newmat(int i, int j ) // return allocated matrix  /* FUNCTION */
 {
 matrix x;
 x.nrows = i;
 x.ncols = j;
 x.data = new double[ i*j ];
 for ( int k = 0 ; k < i*j ; k++ )
 	x.data[k] = 0.0;
 return x;
 }
 
 matrix apply_elwise(matrix x, double f(double)) 
  {
  matrix ans = x ;
  int nr = ans.rows();
  int nc = ans.cols();
  for (int i = 0 ; i < nr; i++ )
    for (int j = 0 ; j < nc; j++ )
     {
     double tmp = f(x.el(i,j));
     set_el(ans,i,j) = tmp;
     }
  return ans;
  }
 
 matrix dapply_elwise(matrix y, matrix mu, double f(double,double))
  {
 // should be bulletproofed, intended only for QL on 2 1-ds
  matrix ans = y;
  int nr = ans.rows();
  int nc = ans.cols();
  for (int i = 0 ; i < nr; i++ )
    for (int j = 0 ; j < nc; j++ )
     {
     double tmp = f(y.el(i,j),mu.el(i,j));
     set_el(ans,i,j) = tmp;
     }
  return ans;
  }
  
 
 matrix mat11(double x)  /* make a matrix object from a scalar */
 {
 matrix go = newmat(1,1);
 set_el(go,0,0) = x;
 return go;
 }
 
 matrix mat11(int x)
 {
 return mat11((double)x);
 }
 
 //void mprint(matrix x)
 //{
 //int nr = x.rows();
 //int nc = x.cols();
 //for (int i = 0 ; i < nr ; i++ )
 // {
 // for (int j = 0 ; j < nc ; j++ )
 //  {printf("%6.3lf ",x.el(i,j));}
 // printf("\n");
 // }
 //}
 
 //static char SCID[] = "@(#) Msort_on_col.cc %I% %D%";
  
 
           static    int Mcpp_compare(const void *i, const void *j)
           {
 // this is explicit type conversion (see pp52-53 of stroustrup 2e
 // necessary for proper behavior of stdlib.h compliant compare in qsort
 		double* x = (double*) i;
 		double* y = (double*) j;
 //               cout << *x;
                int out = 0;
 		if (*x < *y) out = -1;
 		else if (*x > *y) out = 1;
 		return out;
           }
 
 matrix make_row( double in1, double in2, double in3, double in4)
 	{
 	matrix x = newmat(1,4);
 	set_el(x,0,0) = in1;
 	set_el(x,0,1) = in2;
 	set_el(x,0,2) = in3;
 	set_el(x,0,3) = in4;
 	return x;
 	};
 
 matrix transp( matrix base )  /* note -- a preferred way is to dynamically */ /* FUNCTION */
 				/* modify the operation of el */
 {
 matrix x;
 double *look, *load;
 int i, j;
 
 x = newmat( base.cols() , base.rows() );
 load = x.mathead();
 look = base.mathead();
 for ( i = 0 ; i < x.rows() ; i++ )
 	{
 	for ( j = 0 ; j < x.cols() ; j++ )  /* load row-major */
 		{
 		*(load++) = base.el(j,i);
 		}
 	}
 return x;
 }
 matrix operator++( matrix arg )
 	{
 	return( transp(arg) );
 	}
 
 
 
 matrix seq( int start, int end, int gran=1 )  /* FUNCTION */
 {
 if ( gran <= 0 && (start < end) ) error_signal(mcpp_env, SEQ_FAIL_INFINITE_SEQ);
 if ( gran >= 0 && (start > end) ) error_signal(mcpp_env, SEQ_FAIL_INFINITE_SEQ);
 
 int nel = ((int)(abs(end - start)/abs(gran))) + 1;
 
 matrix tmp = newmat( 1 , nel );
 
 double* head = tmp.mathead();
 
 if ( start == end ) *head = (double) start;
 
 else if ( start < end )
 	{
 	for ( ; start <= end ; *(head++) = (double) start , start+=gran );
 	}
 else	
 	{
 	for ( ; start >= end ; *(head++) = (double) start, start+=gran );
 	}
 
 return tmp;
 }
 
 matrix Msort_on_col( matrix& x , int col )   /* col zerobased */
 {
 matrix null;
 if ( col > x.cols() )
 	{
 #ifndef FOR_S
 	error_signal(mcpp_env, MSORT_FAIL_NO_COL);
 #else
 	return(null);
 #endif
 	}
 matrix sortcol = x.submat(0,mat11(col));
 matrix tmp = sortcol || x;
 int wid = tmp.cols() * 8;
 int nel = x.rows();
 double* base = tmp.mathead();
 // this works because compare function looks only at first double of each
 // record, the place where we stuck the key column 
 qsort( (char *)base , nel , wid , Mcpp_compare );
 // return matrix after key column removed
 tmp = tmp.submat(0,seq(1,tmp.cols()-1,1));
 return tmp;
 }
 
 matrix band( matrix& inmat , int bandw )  /** band any matrix to bandwith bandw **/
 {
 	matrix tmp = inmat;
 	for ( int i = 0 ; i < inmat.rows() ; i++ )
 		{
 		for ( int j = i+bandw ; j < inmat.cols() ; j++ )
 			{
 			set_el( tmp , i , j ) = (double)0.;
 			if ( ( i < inmat.cols() ) && ( j < inmat.rows() ) )
 				{
 				set_el( tmp , j , i ) = (double)0.;
 				}
 			}
 		}
 	return tmp;
 }
 
 matrix chol(matrix x)  /* form cholesky decomp of x */  /* FUNCTION */
 			/* x = u'u with u upper right triangular */
 {
 matrix null;
 int p = x.rows();
 
 if ( p != x.cols() )
 	{
 	error_signal(mcpp_env, CHOL_FAIL_NOT_SQR);
 	}
 
 matrix tmp = x;
 
 //double* tmphead = tmp.mathead();
 
 #define loctmp( i , j ) *(tmphead + ( i * p ) + j )  /* element reference for assn */
 
 for ( int i = 0 ; i < p ; i++ )
 	{
 	double accum = 0.;
 
 	for ( int k = 0 ; k <= i-1 ; k++ )
 		{
 		accum += tmp.el(k,i) * tmp.el(k,i);
 		}
 
 	if ( ( tmp.el(i,i) - accum ) >= 0. )
 		{
 		set_el(tmp, i , i ) = (double) sqrt( tmp.el(i,i) - accum );
 		}
 	else 
 		{
 	error_signal(mcpp_env, CHOL_FAIL_ILLG_DIAG);
 		}
 
 
 	for ( int j = i+1  ; j < p ; j++ )
 		{
 		double accum2 = 0.;
 		for ( int k = 0 ; k <= i-1 ; k++ )
 			{
 			accum2 += tmp.el(k,i) * tmp.el(k,j);
 			}
 		if ( tmp.el(i,i) > 0. ) set_el(tmp, i , j ) = ( tmp.el(i,j) - accum2 )/tmp.el(i,i);
 		else 
 			{
 			error_signal(mcpp_env, CHOL_FAIL_NON_POS);
 			}
 	set_el(tmp, j , i ) = 0.;
 		}
 	}
 return tmp;
 }
 
 /* cluscount -- returns integer length of run length encoding */
 
 int cluscount( matrix& disc )
 {
 if (disc.cols() != 1)
 	{
 	error_signal(mcpp_env, CLUSCNT_FAIL_SPLIT_ERR);
 	}
 
 int k = 0;
 
 #define MEL( a, b, c ) a.el( b, c )
 
 int istart = (int)MEL( disc , 0 , 0 );
 int start = 0;
 int end = 0;
 for ( int i = 1 ; i <= disc.rows() ; i++ )
         {
         if (( MEL( disc , i, 0 ) != istart ) ||
                         i == (disc.rows() ) )
                 {
                 k++;
                 start = end+1;
                 istart = (int)MEL( disc, i, 0 );
                 }
         if (start < disc.rows() ) end++ ;
         }
 /* DOES NOT CLEAN */
 return k;
 }
 
 
 matrix col_1s( int i )  /* FUNCTION */
 {
 matrix tmp ;
 
 tmp = newmat( i , 1 );
 double* head = tmp.mathead();
 
 for ( int j = 0 ; j < i ; j++ ) 
 	*(head++) = 1.;
 return tmp;
 }
 
 
  int is_square( matrix& inmat )  /* FUNCTION */
  	{
  	return ( inmat.rows() == inmat.cols() );
  	}
 
 int length( matrix& in )
 	{
 	return in.rows() * in.cols();
 	}
 
 double max( matrix& in )  /* FUNCTION */
 	{
      double mx = in.el(0,0);
 	for ( int i = 0 ; i < in.rows() ; i++ )
 		{
 		for ( int j = 0 ; j < in.cols() ; j++ )
 			{
 			if (in.el(i,j) > mx ) mx = in.el(i,j);
 			}
 		}
 	return mx;
 	}
 
 int compare( const void* i , const void* j )
 {
 int dir = 1;
 
 if ((double *)i < (double *)j) dir = -1;
 else if ((double *)i == (double *)j) dir = 0;
 return dir;
 }
 
 matrix operator||(matrix arg1, matrix arg2) /* horizontal concatenation */  /* FUNCTION */
 	{
 matrix concat;	
 matrix null;
 double *matlook;
 
 	if ( arg1.rows() == 0 )   /* matrix uninitialized */
 		{
 		concat = arg2;
 		return concat;
 		}
 
 	if  ( arg1.rows() != arg2.rows() )
 		{
 		error_signal(mcpp_env, CONCAT_FAIL_ROW_AGRMNT);
 		}
 
 	concat = newmat( arg1.rows() , arg1.cols()+arg2.cols() );
 	matlook = concat.mathead();
 
 	int a1r = arg1.rows();
 	int a1c = arg1.cols();
 	int a2c = arg2.cols();
 
 	for ( int i = 0; i < a1r ; i++ )
 		{
 		for ( int j = 0; j < a1c + a2c ; j++ )
 			{
 			*(matlook++) = ( j < a1c ) ? arg1.el(i,j) : arg2.el(i,j-a1c);
 			}
 		}
 	return concat;
 	}
 
 matrix operator/(matrix& arg1, matrix& arg2) /* vertical concatenation */  /* FUNCTION */
 {
 matrix null;
 matrix concat;	
 double *matlook;
 
 	if ( arg1.cols() == 0 )
 		{
 		concat = arg2;
 		return concat;
 		}
 
 	if  ( arg1.cols() != arg2.cols() )
 		{
 		error_signal(mcpp_env, CONCAT_FAIL_COL_AGRMNT);
 		}
 
 	concat = newmat( arg1.rows()+arg2.rows() , arg1.cols() );
 	matlook = concat.mathead();
 
 	int a1r = arg1.rows();
 	int a2r = arg2.rows();
 	int a1c = arg1.cols();
 	int a2c = arg2.cols();
 
 	for ( int i = 0; i < a1r+a2r ; i++ )
 		{
 		for ( int j = 0; j < a1c ; j++ )
 			{
 			*(matlook++) = ( i < a1r ) ? arg1.el(i,j) : arg2.el(i-a1r,j);
 			}
 		}
 	return concat;
 }
 
 matrix corner( matrix& inmat , int nr , int nc )
 	{
 	matrix out = newmat( nr , nc );
 	for ( int i = 0 ; i < nr ; i++ )
 		{
 		for ( int j = 0 ; j < nr ; j++ )
 			{
 			set_el( out , i , j ) = inmat.el( i , j );
 			}
 		}
 	return out;
 	}
 
 matrix rowseg( matrix& base , int rowstart , int nrows )  /* FUNCTION */
 {
 matrix null;
 matrix x;
 double *look, *load;
 int i, j;
 
 if ( rowstart + nrows > base.rows() || rowstart > base.rows() )
 	{
 	error_signal(mcpp_env, ROWSEG_FAIL_DIM_AGRMT);
 	}
 
 x = newmat( nrows , base.cols() );
 load = x.mathead();
 look = base.mathead() + (rowstart)*base.cols();
 for ( i = 0 ; i < nrows ; i++ )
 	{
 	for ( j = 0 ; j < base.cols() ; j++ )
 		{
 		*(load++) = *(look++);
 		}
 	}
 return x;
 }
 
 matrix delete_row( matrix X, int i)  /* FUNCTION */
 	/* i is row to delete , zerobased */
 {
 matrix result;
 matrix null;
 int ncop, nresel, skip;
 double *resbase, *Xbase;
 
 if (X.rows() <= 1)
 	{
 	error_signal(mcpp_env, DELROW_FAIL_NO_ROW);
 	}
 
 result = newmat( X.rows()-1, X.cols() );
 
 nresel = result.rows() * result.cols() ;
 
 skip = i*X.cols();
 
 resbase = result.mathead();
 Xbase = X.mathead();
 
 for ( ncop = 0 ; ncop < nresel ; ncop++ )
 	{
 	if ( ncop == skip ) Xbase += X.cols();
 	*resbase++ = *Xbase++ ;
 	}
 
 return result;
 }
 
 matrix delete_col( matrix X, int i)  /* FUNCTION */
 	/* i is col to delete , zerobased */
 {
 matrix result;
 matrix null;
 int ncop, nresel, skip;
 double *resbase, *Xbase;
 
 if (X.cols() <= 1)
 	{
 	error_signal(mcpp_env, DELCOL_FAIL_NO_COL);
 	}
 
 result = transp ( delete_row( transp( X ) , i ) );
 
 return result;
 }
 
 matrix extract_diag( matrix inmat )   /* FUNCTION */
 {  /* note, converts to column matrix */
    /* apply form_diag to get actual diag mat */
 matrix null;
 if ( !is_square(inmat) )
 	{
 	error_signal(mcpp_env, EXTRACT_FAIL_NOT_SQR);
 	}
 	
 matrix tmp = newmat( inmat.rows(), 1 );
 double* tmphead = tmp.mathead();
 
 for ( int i = 0 ; i < inmat.rows() ; i++ )
 	{
 	*(tmphead++) = inmat.el(i,i);
 	}
 return tmp;
 }
 
 matrix form_diag( matrix invec )
 {
 matrix null;
 if ( (invec.rows() > 1) && (invec.cols() > 1) )
 	{
 	error_signal(mcpp_env, FORM_FAIL_NOT_VECT);
 	}
 
 int dim = invec.rows() * invec.cols();
 matrix curr = newmat( dim, dim );
 
 for ( int i = 0 ; i < dim ; i++ )
 	{
 	set_el( curr , i , i ) = invec.el(0,i);
 	}
 return curr;
 }
 
 matrix ident( int i )  /* FUNCTION */
 {
 matrix tmp;
 
 tmp = newmat( i , i );
 double* head = tmp.mathead();
 
 for ( int j = 0 ; j < i ; j++ )
 	*(head+( (i+1)*j )) = 1.;
 return tmp;
 }
 
 void isweep( matrix& mat ) /* FUNCTION */
 
 /* algorithm follows Goodnight, Amer Stat Aug 1979 p.149 */
 /* isweep = IN PLACE (avoid copy) */
 
 {
 
 int k , j , i ;
 double d , b, *out, *curtmp;
 
 int tcol = mat.cols();
 int trow = mat.rows();
 
 out = mat.mathead();
 
 for ( k = 0 ; k < trow ; k++ )   /* do 3 */
 
 	{
 
 	d = mat.el(k,k);
 
 	for ( j = 0 ; j < trow ; j++ )   /*do 1 */
 
 		{
 
 		curtmp = out + j + k*(tcol);  /* address to load */
 		*curtmp = mat.el(k,j) / d ;
 
 		}   /* 1 */
 
 	for ( i = 0 ; i < trow ; i++ )  /* do 2 */
 
 		{
 
 		if ( i != k ) 
 
 			{
 
 			b = mat.el(i,k);
 
 			for ( j = 0 ; j < trow ; j++ )
 
 				{
 
 				curtmp = out + j + i*(tcol);
 				*curtmp = mat.el(i,j) - b*mat.el(k,j); 
 
 				}
 
 			curtmp = out + k + i*(tcol);
 			*curtmp = -b / d;
 
 			}   /* end else: 2 */  
 
 		}					/* 2 */	
 
 		curtmp = out + k + k*(tcol);
 		*curtmp = (double) 1.0 / d ;
 
 	}  /* 3 */
 
 }
 
 void kill( matrix& mat )
 	{
 	delete mat.data;
 	mat.ncols = 0;
 	mat.nrows = 0;
 	}
 
 matrix make_row( double in)
 	{
 	matrix x = newmat(1,1);
 	set_el(x,0,0) = in;
 	return x;
 	};
 matrix make_row( double in1, double in2)
 	{
 	matrix x = newmat(1,2);
 	set_el(x,0,0) = in1;
 	set_el(x,0,1) = in2;
 	return x;
 	};
 matrix make_row( double in1, double in2, double in3)
 	{
 	matrix x = newmat(1,3);
 	set_el(x,0,0) = in1;
 	set_el(x,0,1) = in2;
 	set_el(x,0,2) = in3;
 	return x;
 	};
 
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5)
 	{
 	matrix x = newmat(1,5);
 	set_el(x,0,0) = in1;
 	set_el(x,0,1) = in2;
 	set_el(x,0,2) = in3;
 	set_el(x,0,3) = in4;
 	set_el(x,0,4) = in5;
 	return x;
 	};
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6)
 	{
 	matrix x = newmat(1,6);
 	set_el(x,0,0) = in1;
 	set_el(x,0,1) = in2;
 	set_el(x,0,2) = in3;
 	set_el(x,0,3) = in4;
 	set_el(x,0,4) = in5;
 	set_el(x,0,5) = in6;
 	return x;
 	};
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6, double in7)
 	{
 	matrix x = newmat(1,7);
 	set_el(x,0,0) = in1;
 	set_el(x,0,1) = in2;
 	set_el(x,0,2) = in3;
 	set_el(x,0,3) = in4;
 	set_el(x,0,4) = in5;
 	set_el(x,0,5) = in6;
 	set_el(x,0,6) = in7;
 	return x;
 	};
 
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6, double in7, double in8)
 	{
 	return make_row( in1, in2, in3, in4, in5, in6 ) || make_row( in7, in8 );
 	}
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6, double in7, double in8, double in9)
 	{
 	return make_row( in1, in2, in3, in4, in5, in6 ) || make_row( in7, in8 , in9);
 	}
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6, double in7, double in8, double in9, double in10)
 	{
 	return make_row( in1, in2, in3, in4, in5, in6 ) || make_row( in7, in8 , in9,
 	in10);
 	}
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6, double in7, double in8, double in9, double in10, double in11)
 	{
 	return make_row( in1, in2, in3, in4, in5, in6 ) || make_row( in7, in8 , in9,
 	in10, in11);
 	}
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6, double in7, double in8, double in9, double in10, double in11,
 double in12)
 	{
 	return make_row( in1, in2, in3, in4, in5, in6 ) || make_row( in7, in8 , in9,
 	in10, in11, in12);
 	}
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6, double in7, double in8, double in9, double in10, double in11,
 double in12, double in13)
 	{
 	return make_row( in1, in2, in3, in4, in5, in6 ) || make_row( in7, in8 , in9,
 	in10, in11, in12) || make_row( in13 );
 	}
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6, double in7, double in8, double in9, double in10, double in11,
 double in12, double in13, double in14)
 	{
 	return make_row( in1, in2, in3, in4, in5, in6 ) || make_row( in7, in8 , in9,
 	in10, in11, in12) || make_row( in13, in14 );
 	}
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6, double in7, double in8, double in9, double in10, double in11,
 double in12, double in13, double in14, double in15)
 	{
 	return make_row( in1, in2, in3, in4, in5, in6 ) || make_row( in7, in8 , in9,
 	in10, in11, in12) || make_row( in13, in14 , in15);
 	}
 matrix make_row( double in1, double in2, double in3, double in4, double
 in5, double in6, double in7, double in8, double in9, double in10, double in11,
 double in12, double in13, double in14, double in15, double in16)
 	{
 	return make_row( in1, in2, in3, in4, in5, in6 ) || make_row( in7, in8 , in9,
 	in10, in11, in12) || make_row( in13, in14 , in15, in16);
 	}
 
 matrix matlog( matrix& mat )  /* someday, log will just be overloaded.  not now */  /* FUNCTION */
 {
 matrix null;
 double *load;
 matrix tmp;
 
 tmp = newmat( mat.rows() , mat.cols() );
 
 load = tmp.mathead();
 
 for ( int i = 0 ; i < mat.rows() * mat.cols() ; i++ )
 	{
 	if ( mat.data[i] > 0. )
 		*(load++) = log( mat.data[i] );
 	else
 		{
 		error_signal(mcpp_env, MATMATH_FAIL_NEG_MATRIX);
 		}
 	}
 	return tmp;
 }
 
 matrix matexp( matrix& mat )  /* FUNCTION */
 {
 double *load;
 matrix tmp;
 
 tmp = newmat( mat.rows() , mat.cols() );
 
 load = tmp.mathead();
 
 for ( int i = 0 ; i < mat.rows() * mat.cols() ; i++ )
 	{
 	*(load++) = exp( mat.data[i] );
 	}
 return tmp;
 }
 
 matrix matsqrt( matrix& mat )  /* FUNCTION */
 {
 double *load;
 matrix tmp;
 
 tmp = newmat( mat.rows() , mat.cols() );
 
 load = tmp.mathead();
 
 for ( int i = 0 ; i < mat.rows() * mat.cols() ; i++ )
 	{
 	*(load++) = sqrt( mat.data[i] );
 	}
 return tmp;
 }
 
 
 matrix matabs( matrix& mat )  /* FUNCTION */
 {
 double *load;
 matrix tmp;
 
 tmp = newmat( mat.rows() , mat.cols() );
 
 load = tmp.mathead();
 
 for ( int i = 0 ; i < mat.rows() * mat.cols() ; i++ )
 	{
 	*(load++) = fabs( mat.data[i] );
 	}
 return tmp;
 }
 
 
 double min( matrix& in )  /* FUNCTION */
 	{
 	double min = in.el(0,0);
 	for ( int i = 0 ; i < in.rows() ; i++ )
 		{
 		for ( int j = 0 ; j < in.cols() ; j++ )
 			{
 			if (in.el(i,j) < min ) min = in.el(i,j);
 			}
 		}
 	return min;
 	}
 
 double elsum( matrix mat )  /* FUNCTION */
 {
 double sum = 0.;
 for ( int i = 0 ; i < mat.rows() * mat.cols() ; i++ )
 	{
 	sum +=  mat.data[i] ;
 	}
 return sum;
 }
 
 matrix sweep( matrix mat ) /* FUNCTION */
 
 /* algorithm follows Goodnight, Amer Stat Aug 1979 p.149 */
 
 {
 
 int k , j , i ;
 double d , b, *in, *out, *curtmp;
 matrix temp;
 
 temp = newmat(mat.rows(),mat.cols());
 int tcol = temp.cols();
 int trow = temp.rows();
 
 in = mat.mathead();
 out = temp.mathead();
 
 for ( i = 0 ; i < trow*tcol ; /* copy to out */
 		*(out++) = *(in++), i++ )
 		;
 
 out = temp.mathead();
 
 for ( k = 0 ; k < trow ; k++ )   /* do 3 */
 
 	{
 
 	d = temp.el(k,k);
 
 	for ( j = 0 ; j < trow ; j++ )   /*do 1 */
 
 		{
 
 		curtmp = out + j + k*(tcol);  /* address to load */
 		*curtmp = temp.el(k,j) / d ;
 
 		}   /* 1 */
 
 	for ( i = 0 ; i < trow ; i++ )  /* do 2 */
 
 		{
 
 		if ( i != k ) 
 
 			{
 
 			b = temp.el(i,k);
 
 			for ( j = 0 ; j < trow ; j++ )
 
 				{
 
 				curtmp = out + j + i*(tcol);
 				*curtmp = temp.el(i,j) - b*temp.el(k,j); 
 
 				}
 
 			curtmp = out + k + i*(tcol);
 			*curtmp = -b / d;
 
 			}   /* end else: 2 */  
 
 		}					/* 2 */	
 
 		curtmp = out + k + k*(tcol);
 		*curtmp = (double) 1.0 / d ;
 
 	}  /* 3 */
 
 return temp;
 
 }
 
 matrix operator^(matrix& arg1, int arg2)  /* FUNCTION */
 	{
 matrix answer;	
 
 	if  ( arg1.rows() != arg1.rows() )
 		{
 		error_signal(mcpp_env, MATPOW_FAIL_NOT_SQR);
 		}
 
 	if ( arg2 == -1 ) return(sweep( arg1 ));
 	else
 		{
 		error_signal(mcpp_env, MATPOW_FAIL_INVERSE_ONLY);
 		}
 	}
 
 matrix operator-(matrix& arg1, matrix& arg2)  /* FUNCTION */
 	{
 matrix diff;	
 double *matlook;
 
 	if ( ( arg1.rows() != arg2.rows() ) ||
 	     ( arg1.cols() != arg2.cols() ) )
 		{
 		error_signal(mcpp_env, MINUS_FAIL_DIM_AGRMNT);
 		}
 
 	diff = newmat( arg1.rows() , arg1.cols() );
 	matlook = diff.mathead();
 
 	for ( int i = 0; i < arg1.rows() ; i++ )
 		{
 		for ( int j = 0; j < arg1.cols() ; j++ )
 			{
 			*(matlook++) = arg1.el(i,j) - arg2.el(i,j);
 			}
 		}
 	return diff;
 	}
 
 matrix operator-(double arg1, matrix& arg2)  /* FUNCTION */
 	{
 	matrix sum = newmat( arg2.rows(), arg2.cols() );
 
 	double* matlook = sum.mathead();
 
 	for ( int i = 0; i < arg2.rows() ; i++ )
 		{
 		for ( int j = 0; j < arg2.cols() ; j++ )
 			{
 			*(matlook++) = arg1 - arg2.el(i,j);
 			}
 		}
 	return sum;
 	}
 
 void operator-=(matrix& arg1, matrix& arg2)  /* FUNCTION */
 	{
 
 double *matlook;
 
 	if ( ( arg1.nrows != arg2.nrows ) ||
 	     ( arg1.ncols != arg2.ncols ) )
 		{
 #ifndef FOR_S
 		cerr << " matrix-= encounters fatal error:\n";
 		cerr << " matrix-=:  arg1 is ( " << arg1.nrows <<
 				" x " << arg1.ncols << " ) \n";
 		cerr << " matrix-=:  arg2 is ( " << arg2.nrows <<
 				" x " << arg2.ncols << " ) \n";
 		cerr << " matrix-=:  arg1 and arg2 must agree in dimensionality.\n";
 		error_signal(mcpp_env, MINUSEQ_FAIL_DIM_AGRMNT);
 #else
 		fprintf( stderr, "MC++: -= encounters improp dimensioned args. Die.\n");
 		return;
 #endif
 		}
 
 	matlook = arg1.mathead();
 
 	for ( int i = 0; i < arg1.rows() ; i++ )
 		{
 		for ( int j = 0; j < arg1.cols() ; j++ )
 			{
 			*(matlook++) -= arg2.el(i,j) ;
 			}
 		}
 	/* return arg1;  unnecessary */
 	}
 
 matrix operator*(matrix arg1, matrix arg2)  /* FUNCTION */
 {
 matrix prod; 
 double *matlook;
 
 	if  ( arg1.cols() != arg2.rows() )
 		{
 		error_signal(mcpp_env, MULT_FAIL_DIM_AGRMNT);
 		}
 
 	prod = newmat( arg1.rows() , arg2.cols() );
 
 	matlook = prod.mathead();
 	double* a1head = arg1.mathead();
 	double* a2head = arg2.mathead();
 	double* a1cur = a1head;
 	double* a2cur = a2head;
 
 	int a2c = arg2.cols();
 	int a1c = arg1.cols();
 	double* a2h = arg2.mathead();
 
 	for ( int j = 0; j < prod.rows() ; j++ )
 		{
 		for ( int i = 0; i < prod.cols() ; i++ )
 			{
 			a1cur = a1head;
 			a2cur = a2head;
 			for ( int k = 0 ; k < arg2.rows() ; k++ )
 				{
 				*(matlook) += *(a1cur++) * *a2cur;
 				a2cur += a2c;
 				}
 			matlook++;
 			a2head++;
 			}
 		a1head += a1c;
 		a2head = a2h;
 		}
 	return prod;
 }
 
 matrix operator*(double arg1, matrix arg2)  /* FUNCTION */
 {
 matrix prod;
 double *matlook;
 
 prod = newmat( arg2.rows(), arg2.cols() );
 matlook = prod.mathead();
 
 for ( int i = 0 ; i < arg2.rows()*arg2.cols() ; i++ )
 	{
 	*(matlook+i) = arg1 * arg2.data[i] ;
 	}
 return prod;
 }
 
 matrix operator*(matrix arg1, double arg2)  /* FUNCTION */
 {
 return arg2*arg1;
 }
 
 void plug( matrix plugm, matrix socket, int row, int col)  /* FUNCTION */
 {
 int pcol = plugm.cols();
 int prow = plugm.rows();
 
 if (( pcol+col > ( socket.cols() ) ) || ( prow + row > socket.rows() ))
 	{
 	error_signal(mcpp_env, PLUG_FAIL_DIM_AGRMNT);
 	}
 
 double* sockload = socket.mathead() + col + row*(socket.cols());
 double* plughead = plugm.mathead();
 double* sockrow_start = sockload;
 
 for ( int i = 0 ; i < prow ; i++ )
 	{
 	sockload = sockrow_start;
 	for ( int j = 0 ; j < pcol ; j++ )
 		{
 		*(sockload++) = *(plughead++);
 		}
 	sockrow_start += socket.cols();
 	}
 }
 
 matrix operator+(matrix arg1, matrix arg2)  /* FUNCTION */
 	{
 matrix sum;	
 double *matlook;
 
 	if ( ( arg1.rows() != arg2.rows() ) ||
 	     ( arg1.cols() != arg2.cols() ) )
 		{
 		error_signal(mcpp_env, PLUS_FAIL_DIM_AGRMNT);
 		}
 
 	sum = newmat( arg1.rows() , arg1.cols() );
 	matlook = sum.mathead();
 
 	for ( int i = 0; i < arg1.rows() ; i++ )
 		{
 		for ( int j = 0; j < arg1.cols() ; j++ )
 			{
 			*(matlook++) = arg1.el(i,j) + arg2.el(i,j);
 			}
 		}
 	return sum;
 	}
 
 matrix operator+(double arg1, matrix& arg2)  /* FUNCTION */
 	{
 	matrix sum = newmat( arg2.rows(), arg2.cols() );
 
 	double* matlook = sum.mathead();
 
 	for ( int i = 0; i < arg2.rows() ; i++ )
 		{
 		for ( int j = 0; j < arg2.cols() ; j++ )
 			{
 			*(matlook++) = arg1 + arg2.el(i,j);
 			}
 		}
 	return sum;
 	}
 
 void operator+=(matrix& arg1, matrix& arg2)  /* FUNCTION */
 	{
 
 double *matlook;
 	if ( arg1.rows() == 0 ) 
 		{
 		arg1 = arg2;
 		return;
 		}
 
 	if ( ( arg1.rows() != arg2.rows() ) ||
 	     ( arg1.cols() != arg2.cols() ) )
 		{
 		error_signal(mcpp_env, PLUSEQ_FAIL_DIM_AGRMNT);
 		}
 
 	matlook = arg1.mathead();
 
 	for ( int i = 0; i < arg1.rows() ; i++ )
 		{
 		for ( int j = 0; j < arg1.cols() ; j++ )
 			{
 			*(matlook++) += arg2.el(i,j) ;
 			}
 		}
 	/* return arg1;  unnecessary */
 	}
 
 int split( matrix& mat , matrix& disc , matrix recv[] )    /* FUNCTION */
 
 {
 int i=0 , j=0 , k=0, start=0, end=0, len=1;
 double istart;
 
 if (disc.cols() != 1)
 	{
 	error_signal(mcpp_env, SPLIT_FAIL_NOT_COLVEC);
 	}
 
 if (disc.rows() != mat.rows())
 	{
 	error_signal(mcpp_env, SPLIT_FAIL_DIM_AGRMNT);
 	}
 
 k = 0 ;
 
 istart = disc.el(0,0);
 
 for ( i = 1 ; i <= disc.rows() ; i++ )
 	{
 	if ((  disc.el( i , 0 ) != istart ) || i == (disc.rows()) )
 		{
 		len = end - start + 1;
 		recv[k] = rowseg( mat , start , len );
 		k++;
 		start = end + 1;
 		istart =  disc.el( i , 0 );
 		}
 	if (start < disc.rows()) end++;
 	}
 return k;
 }
 
 matrix* split( matrix& mat , matrix& disc )    /* FUNCTION */
 
 {
 int i=0 , j=0 , k=0, start=0, end=0, len=1;
 double istart;
 matrix *recv;
 
 if (disc.cols() != 1)
 	{
 	error_signal(mcpp_env, SPLIT_FAIL_NOT_COLVEC);
 	}
 
 if (disc.rows() != mat.rows())
 	{
 	error_signal(mcpp_env, SPLIT_FAIL_DIM_AGRMNT);
 	}
 
 k = 0 ;
 
 int nclus = cluscount( disc );
 if (!(recv = (matrix *)calloc( nclus, (unsigned)sizeof(class matrix))))
 	{
 	error_signal(mcpp_env, SPLIT_FAIL_NO_POINTERS);
 	}
 
 istart = disc.el(0,0);
 
 for ( i = 1 ; i <= disc.rows() ; i++ )
 	{
 	if ((  disc.el( i , 0 ) != istart ) || i == (disc.rows()) )
 		{
 		len = end - start + 1;
 		recv[k] = rowseg( mat , start , len );
 		k++;
 		start = end + 1;
 		istart =  disc.el( i , 0 );
 		}
 	if (start < disc.rows()) end++;
 	}
 return recv;
 }
 
 matrix matrix::submat( matrix row_req , int dummy )  /* FUNCTION */
 	{
 		/* set up the (complete) col req */
 	int srccols = this->ncols;
 	matrix crq = newmat( 1 , srccols );
 
 	double* load = crq.mathead();
 
 	for ( int j = 0 ; j < srccols ; j++ )
 		{
 		*(load++) = (double)j;
 		}
 
 	matrix tmp = this->submat( row_req, crq );
 	return tmp;
 	}
 
 matrix matrix::submat( int dummy , matrix col_req )  /* FUNCTION */
 	{
 	int srcrows = this->nrows;
 	matrix rrq = newmat( 1 , srcrows ) ;
 
 	double* load = rrq.mathead();
 
 	for ( int i = 0 ; i < srcrows ; i++ )
 		{
 		*(load++) = (double)i;
 		}
 	matrix tmp = this->submat( rrq, col_req );
 	return tmp;
 	}
 
 matrix matrix::submat( matrix row_req , matrix col_req )  /* FUNCTION */
 	{
 matrix crq, rrq;
 matrix null;
 	/* put "request mats " in canonical form -- as row vecs */
 	if ( row_req.rows() == 1 ) /* is a row vector */
 		{
 		rrq = row_req;
 		}
 	else if ( row_req.cols() == 1 )
 		{
 		rrq = transp( row_req );
 		}
 	else
 		{
 		error_signal(mcpp_env, SUBMAT_FAIL_NOT_ROWVEC);
 		}
 	if ( col_req.rows() == 1 ) /* is a row vector */
 		{
 		crq = col_req;
 		}
 	else if ( col_req.cols() == 1 )
 		{
 		crq = transp( col_req );
 		}
 	else
 		{
 		error_signal(mcpp_env, SUBMAT_FAIL_NOT_COLVEC);
 		}
 	int n_r_pulls = rrq.cols();
 	int n_c_pulls = crq.cols();
 	int srccols = this->ncols;
 	int srcrows = this->nrows;
 
 	if ( (int)max(rrq) > (srcrows-1))
 		{
 		error_signal(mcpp_env, SUBMAT_FAIL_CHECK_ROWS);
 		}
 	if ( (int)max(crq) > (srccols-1))
 		{
 		error_signal(mcpp_env, SUBMAT_FAIL_CHECK_COLS);
 		}
 
 	matrix tmp = newmat( n_r_pulls, n_c_pulls );
 
 	double* load = tmp.mathead();
 
 	for ( int i = 0 ; i < n_r_pulls ; i++ )
 		{
 		for ( int j = 0 ; j < n_c_pulls ; j++ )
 			{
 			*(load++) = this->data[ (int)rrq.el(0,i)*srccols + (int)crq.el(0,j) ];
 			}
 		}
 	return tmp;
 	}
 matrix symdet( matrix& mat ) /* FUNCTION */
 
 /* algorithm follows Goodnight, Amer Stat Aug 1979 p.149 */
 
 {
 
 int k , j , i ;
 double d , b, *in, *out, *curtmp, detsc = 1.;
 matrix temp, detmat;
 
 temp = newmat(mat.rows(),mat.cols());
 detmat = newmat( 1 , 1 );
 
 in = mat.mathead();
 out = temp.mathead();
 
 for ( i = 0 ; i < mat.rows()*mat.cols() ; /* copy to out */
 		*(out++) = *(in++), i++ )
 		;
 
 out = temp.mathead();
 detsc = temp.el(0,0);
 
 for ( k = 0 ; k < mat.rows() ; k++ )   /* do 3 */
 
 	{
 
 	d = temp.el(k,k);
 
 	for ( j = 0 ; j < mat.rows() ; j++ )   /*do 1 */
 
 		{
 
 		curtmp = out + j + k*(temp.cols());  /* address to load */
 		*curtmp = temp.el(k,j) / d ;
 
 		}   /* 1 */
 
 	for ( i = 0 ; i < mat.rows() ; i++ )  /* do 2 */
 
 		{
 
 		if ( i != k ) 
 
 			{
 
 			b = temp.el(i,k);
 
 			for ( j = 0 ; j < mat.rows() ; j++ )
 
 				{
 
 				curtmp = out + j + i*(temp.cols());
 				*curtmp = temp.el(i,j) - b*temp.el(k,j); 
 
 				}
 
 			curtmp = out + k + i*(temp.cols());
 			*curtmp = -b / d;
 
 			}   /* end else: 2 */  
 
 		}					/* 2 */	
 
 		curtmp = out + k + k*(temp.cols());
 
 			
 		if ( k <  ( mat.rows() -1 ) )
 			detsc = detsc * temp.el( k+1 , k+1 ) ;
 
 		*curtmp = (double) 1.0 / d ;
 
 	}  /* 3 */
 
 curtmp = detmat.mathead();
 *curtmp = detsc;
 return detmat;
 
 }
 
 matrix toeplitz( matrix in ) /* FUNCTION */
 {
 matrix toep, tin;
 int n, p;
 int inrows = in.rows();
 int incols = in.cols();
 
 if ( (inrows > incols) ? inrows % incols : incols % inrows ) // !=0
 	{
 	error_signal(mcpp_env, TOEPLITZ_FAIL_DIM_AGRMNT);
 	}
 
 if ( inrows > incols )   // n.p * p
 	{
 	p = incols;
 	n = inrows/p;
 	tin = in;
 	}
 else			// p * n.p
 	{
 	p = inrows;
 	n = incols/p;
 	tin = transp(in);   /* to permit plucking of boxes */
 	}
 
 toep = newmat( n*p, n*p );
 
 for ( int i = 0 ; i < n ; i++ )
 	{
 	matrix tmp = rowseg( tin, i*p , p );  /* boxes plucked as rowsegs */
 
 	if ( i == 0 )
 		{
 		for ( int j = 0 ; j < n ; j++ )
 			{
 			if ( inrows > incols )
 				plug( tmp, toep, j*p, j*p );
 			else
 				plug( transp(tmp), toep, j*p, j*p );
 			}
 		}
 	else
 		{
 		for ( int j = 0 ; j < n-i ; j++ )
 			{
 		/*	if ( inrows > incols )
 				{ */
 				plug( transp( tmp ), toep , j*p , (j+i)*p);
 				plug(  tmp , toep , (j+i)*p , j*p);
 				/* }
 			else
 				{
 				plug(  tmp , toep , j*p , (j+i)*p);
 				plug(  transp(tmp) , toep , (j+i)*p , j*p);
 				} */
 			}
 		}
 	}
 return toep;
 }
 
 double trace( matrix mat )
 	{
 	double tmp = 0.;
 	for ( int i = 0 ; i <  mat.rows() ; i++ )
 		tmp += mat.el(i,i);
 	return tmp;
 	}
 
 matrix vec2diag( matrix& x )
 	{
 	matrix null;
 	if ( x.cols() != 1 && x.rows() != 1 )
 		{
 		fprintf(stderr,"MC++: vec2diag: need arg with 1 row or 1 col\n");
 #ifndef FOR_S
 		error_signal(mcpp_env, VEC2DIAG_FAIL_NEED_DIM1);
 #else
 		return(null);
 #endif
 		}
 	if ( x.rows() == 1 ) x = transp(x);
 
 	matrix out = newmat( x.rows(), x.rows() );
 
 	for ( int i = 0 ; i < x.rows() ; i++ )
 		{
 		set_el( out, i, i) = x.el(i,0);
 		}
 	return out;
 	}
 
 // start of GEE solution utilities
 
 
 double matmaxabs(matrix x) {
  matrix absm = apply_elwise(x,(double(*)(double))fabs);
  return max(absm);
  }
 
 matrix mult_like_S(matrix x1, matrix x2)
  {
 /* really only takes a conforming column matrix and 
   "replicates it" like Splus does for elwise product */
  int nr1 = x1.rows();
  int nc1 = x1.cols();
  int nr2 = x2.rows();
  int nc2 = x2.cols();
  int confto2 = 0;
  int NC = 0;
  if (nr1 != nr2 ) error_signal(mcpp_env, MULT_FAIL_DIM_AGRMNT);
  if (nc1 == 1 & nc2 >= 1) { NC = nc2; confto2 = 1; }
  else if (nc2 == 1 & nc1 >=1 ) {NC = nc1; confto2 = 0; }
  else error_signal(mcpp_env, MULT_FAIL_DIM_AGRMNT);
  matrix ans = newmat(nr1,NC);
  for (int i = 0 ; i < nr1 ; i++ )
    for (int j = 0 ; j < NC ; j++ )
      {
      double tmp = 0.;
      if (confto2) tmp = x2.el(i,j) * x1.el(i,0);
      else tmp = x1.el(i,j) * x2.el(i,0);
      set_el(ans,i,j) = tmp;
      }
  return ans;
  }
 
 
 matrix from_S( double* x , int xrow, int xcol )  /* FUNCTION */
 {
 matrix tmp = newmat( xcol, xrow );
 double *head = tmp.mathead();
 
 for ( int j = 0 ; j < xrow * xcol ; j++ )
 	{
 	*(head++) = *(x++);
 	}
 tmp = transp(tmp);
 return tmp;
 }
 
 void to_S( matrix inmat , double* target ) /* FUNCTION */
 {
 double *dummy;
 
 dummy = target;
 
 for ( int i = 0 ; i < inmat.cols() ; i++ )
 	{
 	for ( int j = 0 ; j < inmat.rows() ; j++ )
 		{
 		*(dummy++) = inmat.el(j,i);
 		}
 	}
 }
 
 double alogit(double x) {
  double tmp = exp(x);
  return tmp/(1.0+tmp);
  }
 double alogitxcomp(double x) {
  double tmp = alogit(x);
  return tmp*(1-tmp);
  }
 double muxcomp(double x) {
  return x*(1-x);
  }
 double square(double x) {
  return x*x;
  }
 double reciproot(double x) {
  return 1./sqrt(x);
  }
 double recip(double x) {
  return 1./x;
  }
 double qrecip(double x) {
  return 1./(x*x);
  }
 double nqrecip(double x) {
  return -1./(x*x);
  }
 double const1(double x) {
  return 1.;
  }
 double dident(double x) {
  return x;
  }
 double locexp(double x) {
  return exp(x);
  }
 
 // functions of y,mu
 // evaluated 
 // by matrix dapply_elwise(matrix y, matrix mu, double f(double,double))
 // ql_constv
 // ql_identv
 // ql_quadrv
 // ql_bincomplv
 
 double ql_constv(double y, double mu) {
  return -(y - mu)*(y - mu)/2.;
  }
 double ql_identv(double y, double mu) {
  return y*log(mu) - mu;
  }
 double ql_quadrv(double y, double mu) {
  return y/mu - log(mu);
  }
 double ql_bincomplv(double y, double mu) {
  return y*log(mu/(1.-mu)) + log(1.-mu);
  }
 
 //void print_mat(matrix x)
 // {
 // int nr = x.rows();
 // int nc = x.cols();
 // for (int i = 0 ; i < nr; i++ )
 //  {
 //  for (int j = 0 ; j < nc ; j++ )
 //   {
 //   printf("%lf ",x.el(i,j));
 //   }
 //  printf("\n");
 //  }
 // }
 
 matrix zeralp(matrix PRin, matrix ID, matrix TIMin, double phi, int p, matrix alpin,
 double atol, int amaxit)
 {return mat11(0.0);}
 
#ifdef EXPERIMENTAL_USERCOR
 extern matrix alpfun_user(matrix , matrix , matrix , double , int , matrix ,
 double , int );
#endif
 
 matrix LZ_exchalp(matrix PRin, matrix ID, matrix TIMin, double phi, int p, matrix alpin,
 double atol, int amaxit)
 // TIMin, alpin, atol, amaxit included for fixed signature
  {
  int I = cluscount(ID);
  matrix *PR;
  PR = split( PRin, ID );
  double A0 = 0.;
  double den = 0.;
  for (int i = 0 ; i < I ; i++ )
   {
   int ni = PR[i].rows() ;
   if ( ni > 1 )
      {
      A0 = A0 + (elsum(PR[i] * transp(PR[i])) - (transp(PR[i])*PR[i]).el(0,0))/2.;
      den += .5*ni*(ni-1);
      }
   }
  double ans = A0/(phi*(den-(double)p));
  return mat11(ans);
  }
 
 matrix Heidel_unstr(matrix PRin, matrix ID, matrix TIMin, double phi, int p, matrix alpin,
 double atol, int amaxit)
 // TIMin, alpin, atol, amaxit included for fixed signature
 // returns full qxq working covariance est
  {
  int I = cluscount(ID);
  matrix *PR;
  matrix *Ti;
  PR = split( PRin, ID );
  Ti = split( TIMin, ID );
  int maxtime = 0;
  int i = 0;
  for (i = 0 ; i < I; i++ )
   {
   int ti = Ti[i].rows();
   if (ti > maxtime) maxtime = (int) Ti[i].el(ti-1,0);
   }
  int q = maxtime+1; // even if no-one came at all times
  matrix alp = newmat(q,q);
  matrix den = newmat(q,q);
  for (i = 0 ; i < I ; i++ )
   {
   int ni = PR[i].rows() ;
   if ( ni > 1 )
     {
     for (int j = 0 ; j < ni; j++ )
      {
      for (int k = j ; k < ni ; k++ )
        {
        int tj = (int)Ti[i].el(j,0);
        int tk = (int)Ti[i].el(k,0);
        set_el(alp,tj,tk) = alp.el(tj,tk) + PR[i].el(j,0)*PR[i].el(k,0);
        set_el(den,tj,tk) = den.el(tj,tk) + 1.0;
        }
      }
     } 
    } /* done with sums */
   den = den+transp(den);
   alp = alp+transp(alp);
   for (i = 0 ; i < q ; i++ )
     {
     set_el(den,i,i) = (double)I;
     set_el(alp,i,i) = alp.el(i,i)/ 2.;
     }
   for (i = 0 ; i < q ; i++ )
    for (int j = 0 ; j < q ; j++ )
     set_el(alp,i,j) = alp.el(i,j)/ (den.el(i,j)-p);
  return alp;
  }
 
 
 matrix mident(matrix alp, int ni, matrix tim) {
 // alp, tim included for fixed signature
  return ident(ni);
  }
 
#ifdef EXPERIMENTAL_USERCOR
 extern matrix wcorinv_user(matrix , int , matrix );
#endif
 
 matrix exinv(matrix alp, int ni, matrix tim) {
 // tim included just for signature
  double Alp = alp.el(0,0);
  matrix ans = newmat(ni,ni);
  double n = (double)ni;
  double den = (n-1.)*Alp*Alp - (n-2.)*Alp - 1.;
  for (int i = 0 ; i < ni ; i++ )
   {
   set_el(ans,i,i) = -((n-2.)*Alp + 1.)/den;
   for (int j = i+1 ; j < ni ; j++ )
    {
    double tmp = Alp/den;
    set_el(ans,i,j) = tmp;
    set_el(ans,j,i) = tmp;
    }
   }
  return ans;
  }
 
 matrix unstrinv(matrix alp, int ni, matrix tim) {
  matrix ans = newmat(ni,ni);
  for (int i = 0 ; i < ni ; i++ )
   {
   for (int j = 0 ; j < ni ; j++ )
    {
    if (i!=j) set_el(ans,i,j) = alp.el((int)tim.el(i,0),(int)tim.el(j,0)) /
          sqrt(alp.el((int)tim.el(i,0),(int)tim.el(i,0))*
               alp.el((int)tim.el(j,0),(int)tim.el(j,0)));
    }
   }
  for (int i = 0 ; i < ni ; i++ )
        set_el(ans,i,i) = 1.0;
  return sweep(ans);
  }
 
 matrix fominv(matrix alp, int ni, matrix tim) {
 // yes there are analytic forms
 // but not useful for non-lattice val'd obstimes
  matrix ans = newmat(ni,ni);
  double rho = alp.el(0,0);
  for (int i = 0 ; i < ni ; i++ )
   {
   set_el(ans,i,i) = 1.;
   for (int j = i+1 ; j < ni ; j++ )
    {
    double tmp = pow(rho,fabs(tim.el(i,0)-tim.el(j,0)));
    set_el(ans,i,j) = tmp;
    set_el(ans,j,i) = tmp;
    }
   }
  return sweep(ans);
  }
 
 
 double LZ_scalefun( matrix PRin, matrix ID, matrix TIMin, int p )
  {
  double N = (double) PRin.rows();
  return elsum(apply_elwise(PRin,(double(*)(double))square))/(N-(double)p);
  }
  
 
 matrix LDLt(matrix A)
 {
 int n = A.rows();
 matrix v = newmat(n-1,1);
 double s = 0.0;
 for (int k = 0; k < n ; k++ ) /* 10 */
  {
  for (int p = 0; p < k; p++ )
    {   /* 20 */
    double tmp1 = A.el(p,p)*A.el(k,p);
    set_el(v,p,0) = tmp1;
    }  /* end 20 */
  s = 0.0;
  for (int p = 0; p < k; p++ )
    {  /* 25 */
    s = s+A.el(k,p)*v.el(p,0);
    }  /* end 25 */
  double tmp4 = A.el(k,k) - s;
  set_el(A,k,k) = tmp4;
  for (int i = k+1; i < n; i++ )
    {  /* 30 */
    s = 0.0; for (int p = 0; p < k ; p++ )
       {
       s = s + A.el(i,p)*v.el(p,0);
       }
    double tmp2 = (A.el(i,k) - s)/A.el(k,k);
    set_el(A,i,k) = tmp2;
    }  /* end 30 */
  }  /* end 10 */
  return A;
 }
 
 matrix Kmat(int rank)
 {
 /* antidiagonal */
 matrix K = newmat(rank,rank);
 for (int i = 0; i < rank ; i++ )
  set_el(K,i,rank-i-1) = 1.0;
 return K;
 }
 
 matrix lowutri(matrix A)
 {
 int n = A.rows();
 matrix B = newmat(n,n);
 for (int i = 0 ; i < n ; i++ )
  {
  set_el(B,i,i) = 1.0;
  for (int j = i+1; j < n ; j++ )
   set_el(B,j,i) = A.el(j,i);
  }
 return B;
 }
 
 matrix diffL( matrix wi(matrix,int,matrix), matrix alp0, int n, 
          matrix tim, double eps ) 
  {
  int q = alp0.rows();
  matrix alp1 = newmat(q,1);
  matrix R0 = wi( alp0, n, tim );
 
  for (int j = 0 ; j < q ; j++ )
     set_el(alp1,j,0) = alp0.el(j,0) + eps;
 
  matrix R1 = wi( alp1, n, tim );
 
  matrix LR0 = lowutri(LDLt(R0));
  matrix LR1 = lowutri(LDLt(R1));
 
  return (LR1 - LR0)*(1./eps);
  }
 
 matrix diffLK( matrix wi(matrix,int,matrix), matrix alp0, int n, 
          matrix tim, double eps ) 
  {
  int q = alp0.rows();
  matrix alp1 = newmat(q,1);
  matrix K = Kmat(n);
  matrix R0 = K*wi( alp0, n, tim )*K;
 
  for (int j = 0 ; j < q ; j++ )
     set_el(alp1,j,0) = alp0.el(j,0) + eps;
 
 
  matrix R1 = K*wi( alp1, n, tim )*K;
 
  matrix LR0 = K*lowutri(LDLt(R0))*K;
  matrix LR1 = K*lowutri(LDLt(R1))*K;
 
  return (LR1 - LR0)*(1./eps);
  }
 
 
 // diffL( fominv, alp0, 4, dumtim, 0.000001 )
 
 matrix UG_fom_eval(matrix PRin, matrix ID, matrix TIMin, double phi, 
     int p, matrix alpin)
   {
   double rho = alpin.el(0,0);
   double u = 0.;
   double UUt = 0;
   double dUdr = 0.;
   int nc = cluscount(ID);
   matrix* e = split(PRin, ID);
   matrix* d = split(TIMin, ID);
  //printf("UJiter= %d\n",iter);
  //printf("rho= %lf\n",rho);
   for ( int i = 0 ; i < nc ; i++ )
    {
    double ui = 0;
    int ni = e[i].rows();
    for (int j = 1; j < ni; j++) // deliberate late start
     {
     double dij = d[i].el(j,0)-d[i].el(j-1,0);
     double r2d = pow(rho, dij);
     double r22d = pow(rho, 2.*dij);
     double r2dm = pow(rho, dij-1.);
     double r22dm = pow(rho, 2.*dij-1.);
     double eij = e[i].el(j,0);
     double eijm = e[i].el(j-1,0);
     double ujcomp = dij*r2dm*((r2d*eijm - eij)*eijm + (r2d*eij - eijm)*eij)/(1.-r22d);
     double uGcomp = dij*r22dm*((eij-r2d*eijm)*(eij-r2d*eijm) + 
        (eijm-r2d*eij)*(eijm-r2d*eij) - 2.*phi*(1.-r22d))/((1.-r22d)*(1.-r22d));
     u = u + ujcomp + uGcomp;
     ui = ui + ujcomp + uGcomp;
     dUdr += dij*dij*pow(rho,2.*dij-2.)/(1.-r22d);
     }
    UUt += ui*ui;
    }
  dUdr = 2.*phi*dUdr;
  return transp(make_row(u,0.,0.));
  }
 
 matrix UQ_fom_eval(matrix PRin, matrix ID, matrix TIMin, double phi,
     int p, matrix alpin)
   {
   double rho = alpin.el(0,0);
   double u = 0.;
   double UUt = 0;
   double dUdr = 0.;
   int nc = cluscount(ID);
   matrix* e = split(PRin, ID);
   matrix* d = split(TIMin, ID);
  //printf("UJiter= %d\n",iter);
  //printf("rho= %lf\n",rho);
   for ( int i = 0 ; i < nc ; i++ )
    {
    double ui = 0;
    int ni = e[i].rows();
    for (int j = 1; j < ni; j++) // deliberate late start
     {
     double dij = d[i].el(j,0)-d[i].el(j-1,0);
     double r2d = pow(rho, dij);
     double r22d = pow(rho, 2.*dij);
     double r2dm = pow(rho, dij-1.);
     double r22dm = pow(rho, 2.*dij-1.);
     double eij = e[i].el(j,0);
     double eijm = e[i].el(j-1,0);
     double ujcomp = dij*r2dm*((r2d*eijm - eij)*eijm + (r2d*eij - eijm)*eij)/(1.-r22d);
     double uGAMcomp = dij*r22dm*((eij-r2d*eijm)*(eij-r2d*eijm) +
        (eijm-r2d*eij)*(eijm-r2d*eij) )/((1.-r22d)*(1.-r22d));
     u = u + ujcomp + uGAMcomp;
     ui = ui + ujcomp + uGAMcomp;
     dUdr += dij*dij*pow(rho,2.*dij-2.)/(1.-r22d);
     }
    UUt += ui*ui;
    }
  dUdr = 2.*phi*dUdr;
  return transp(make_row(u,0.,0.));
 }
 
 
 matrix UJ_fom_eval(matrix PRin, matrix ID, matrix TIMin, double phi, 
     int p, matrix alpin)
   {
 // proven to obtain s.e. components
 // returns U, dUdr, UUt (summed)
   double rho = 0.0;
   rho = alpin.el(0,0);
 //printf("rho = %lf\n", rho);
   double u = 0.;
 //printf("u = %lf\n", u);
   double dUdr = 0.;
   double UUt = 0;
   int nc = cluscount(ID);
   matrix* e = split(PRin, ID);
   matrix* d = split(TIMin, ID);
   for ( int i = 0 ; i < nc ; i++ )
    {
    double ui = 0;
    int ni = e[i].rows();
    for (int j = 1; j < ni; j++) // deliberate late start
     {
     double dij = d[i].el(j,0)-d[i].el(j-1,0);
     double r2d = pow(rho, dij);
     double r22d = pow(rho, 2.*dij);
     double r2dm = pow(rho, dij-1.);
     double eij = e[i].el(j,0);
     double eijm = e[i].el(j-1,0);
     double utmp = dij*r2dm*((r2d*eijm - eij)*eijm + (r2d*eij - eijm)*eij)/(1.-r22d);
     u = u + utmp;
     ui += utmp;
     dUdr += dij*dij*pow(rho,2.*dij-2.)/(1.-r22d);
     }
    UUt += ui*ui;
    }
   dUdr = 2.*phi*dUdr;
 //printf("u = %lf\n", u);
   return transp(make_row(u,dUdr,UUt));
  }
 
 matrix UJ_equi_eval(matrix PRin, matrix ID, matrix TIMin, double phi, 
     int p, matrix alpin)
   {
   void matd(matrix);
   double rho = alpin.el(0,0);
   double u = 0.;
   double dUdr = 0.;
   double UUt = 0;
   int nc = cluscount(ID);
   matrix* e = split(PRin, ID);
   matrix* d = split(TIMin, ID);
   for ( int i = 0 ; i < nc ; i++ )
    {
    double ui = 0;
    int ni = e[i].rows();
    matrix K = Kmat(ni);
    GetUjStruc( dB, JJ, BB, exinv, alpin, ni, d[i], .000001 )
    GetUjKStruc( dBK, JJK, BBK, exinv, alpin, ni, d[i], .000001 )
    u += (transp(e[i])*dB*JJ*transp(BB)*e[i]).el(0,0) +
             (transp(e[i])*dBK*K*JJK*K*transp(K*BBK*K)*e[i]).el(0,0);
 // you can do the averaging of uj1 and uj2 at this stage (this is
 // just uj2
    UUt += ui*ui;
    }
   return transp(make_row(u,0.,UUt));
  }
 
 
 /* Here is start of revised engine */
 
 /* YUA = yags U_alpha */
 matrix yua_generic_eval( matrix f(matrix,matrix,matrix,double,int,matrix), 
         matrix PRin, matrix ID, matrix TIMin,
         double phi, int p, matrix alp )
  {
  //matrix (*gf)(matrix,matrix,matrix,double,int,matrix);
 //matrix UG_fom_eval(matrix PRin, matrix ID, matrix TIMin, double phi, 
  //   int p, matrix alpin)
  //gf = &f;
  //f = &UG_fom_eval;
  matrix val = (matrix)((*f)(PRin, ID, TIMin, phi, p, alp));
  return val;
  }
 
 #define YUA_FSIG matrix f(matrix,matrix,matrix,double,int,matrix)
 #define YUA_PARMS \
     matrix PRin, matrix ID, matrix TIMin, double phi, int p, matrix alp
 #define YUA_SIG matrix,matrix,matrix,double,int,matrix
 #define YUA_ARGS \
     PRin, ID, TIMin, phi, p, alp
 #define YUA_ARGS1 \
     PRin, ID, TIMin, phi, p, alp1
 #define YUA_ARGS2 \
     PRin, ID, TIMin, phi, p, alp2
 
 matrix yua_grid_minabs( YUA_FSIG, YUA_SIG, double , double , int );
 matrix yua_secant_solve( YUA_FSIG, YUA_SIG, double , double , int , double,
            double, int);
 
 matrix yua_grid_minabs( YUA_FSIG, YUA_PARMS, double lb, double ub, int npts )
  {
  double gran = (ub - lb)/(double)npts;
  double val = 100000.0;
  double ans = 0.;
  double x = lb;
  for (int i = 0 ; i < npts; i++ )
    {
    alp = mat11(x);
    double nval = yua_generic_eval( f, YUA_ARGS ).el(0,0);
    //printf("x=%lf f(x)=%lf\n",x,nval);
    if (fabs(nval)<val) {
       ans = x;
       val = fabs(nval);
       }
    x=x+gran;
    }
  return mat11(ans);
  }
    
 matrix yua_secant_solve( YUA_FSIG, YUA_PARMS, double lb, 
                double ub, int npts, double del , double eps, int maxiter)
  {
 // program assumes that f returns 3-element matrix, only the first
 // of which is relevant to solution.  remaining components
 // are just passed on in addition to the solution 
  matrix tmp;
  double x1 = yua_grid_minabs( f, YUA_ARGS, lb, ub, npts ).el(0,0);
  double x2 = x1 + del;
  double x3 = 0.;
  int iter = 0;
  double f2 = 100.*eps;
  while ((fabs(x1-x2)>eps || fabs(f2) > eps) && iter < maxiter)
   {
   matrix alp1 = mat11(x1);
   matrix alp2 = mat11(x2);
   double f1 = yua_generic_eval( f, YUA_ARGS1 ).el(0,0);
   tmp = yua_generic_eval( f, YUA_ARGS2 );
   f2 = tmp.el(0,0);
 //   printf("x1=%lf f(x1)=%lf\n",x1,f1);
 //  printf("x2=%lf f(x2)=%lf\n",x2,f2);
   x3 = x2 - (f2 * (x2-x1))/(f2-f1);
   x1 = x2;
   x2 = x3;
   }
  return transp(make_row(x2,tmp.el(0,0),tmp.el(1,0),tmp.el(2,0)));
  }
    
 extern "C" {
 
 //jmp_buf mcpp_env;
 
  void yags_engine(int* n, int* p, int* q, double* xin,
 double* yin, double* idin, double* tin, 
 double* b0, double* bout, int* maxiter, 
 double* tol, int* famcode, int* corcode, double* alpin, double* alpout,
 double* phiout, double* bcov_nai, double* bcov_rob, 
 double* Ua, double* dUda, double* sum_uut,
 double* Uatol, int* Uamaxit, int* verbose,
 double* Uagrid_lo, double* Uagrid_hi, int* Uagrid_npts,
 double* Uasecant_del, double* weightin,
 double* QLS, double* AIC_PAN, double* GAU_LL, double* OFFSET,
 int* fixscale)
  {
 /* to match a string "target" in char** strarg, use */
 /* (strcomp(*strarg, (const char *)"target") == 0) */
 
 
 // functions of y,mu
 // evaluated 
 // by matrix dapply_elwise(matrix y, matrix mu, double f(double,double))
 // ql_constv
 // ql_identv
 // ql_quadrv
 // ql_bincomplv
 
 double ql_val, qls_val;
 /* famcodes are 1=gaussian, 2=binom, 3=poisson, 4=gamma, 
      5=log link and quadratic var, 6=id link and lin var,
      7 = id link and quad var, 8=log link and const var, 9=log link and lin variance (sic) */
 const int gaussian = 1;
 const int binomial = 2;
 const int poisson = 3;
 const int gamma = 4;
 const int loglinkquadvar = 5;
 const int idlinklinvar = 6;
 const int idlinkquadvar = 7;
 const int loglinkconstvar = 8;
 const int loglinklinvar = 9;
 
 void* mu;
 void* dmu_deta; // partial mu(eta) / partial eta
 void* variance;
 void* qlfun;
 
 switch (*famcode) {
  printf("famco %d\n", *famcode);
  case gaussian:
    mu = (void *)dident;
    dmu_deta = (void *)const1;
    variance = (void *)const1;
    qlfun = (void *)ql_constv;
    break;
  case binomial:
    mu = (void *)alogit;
    dmu_deta = (void *)alogitxcomp;
    variance = (void *)muxcomp;
    qlfun = (void *)ql_bincomplv;
    break;
  case poisson:
    mu = (void *) locexp; //(void *)exp;
    dmu_deta = (void *) locexp; //(void *)exp;
    variance = (void *)dident;
    qlfun = (void *)ql_identv;
    break;
  case gamma:
    mu = (void *)recip;
    dmu_deta = (void *)nqrecip;
    variance = (void *)square;
    qlfun = (void *)ql_quadrv;
    break;
  case loglinkquadvar:
    mu = (void *) locexp; //(void *)exp;
    dmu_deta = (void *) locexp; //(void *)exp;
    variance = (void *)square;
    qlfun = (void *)ql_quadrv;
    break;
  case idlinklinvar:
    mu = (void *)dident;
    dmu_deta = (void *)const1;
    variance = (void *)dident;
    qlfun = (void *)ql_identv;
    break;
  case idlinkquadvar:
    mu = (void *)dident;
    dmu_deta = (void *)const1;
    variance = (void *)square;
    qlfun = (void *)ql_quadrv;
    break;
  case loglinkconstvar:
    mu = (void *) locexp; //(void *)exp;
    dmu_deta = (void *) locexp; //(void *)exp;
    variance = (void *)const1;
    qlfun = (void *)ql_constv;
    break;
  case loglinklinvar:
    mu = (void *) locexp; //(void *)exp;
    dmu_deta = (void *) locexp; //(void *)exp;
    variance = (void *)dident;
    qlfun = (void *)ql_identv;
    break;
  }  
 
 matrix (*wcorinv)(matrix,int,matrix);
 matrix (*alpfun)(matrix,matrix,matrix,double,int,matrix,double,int); // will get extended signature
 matrix (*alpfun_eval)(YUA_SIG); // 
 
 const int identity = 1;
 const int exchange_LZ = 2;
 const int UJ_fomstr = 3;
 const int UG_fomstr = 4;
 const int UQ_fomstr = 5;
 const int heidel_unstr = 6;
 const int UJ_equistr = 7;
 const int UJ_equimart = 8;
 const int User_cor = 9;
 
 switch (*corcode) {
  case identity:
   alpfun = &zeralp;
   wcorinv = &mident;
   break;
  case exchange_LZ:
   alpfun = &LZ_exchalp;
   wcorinv = &exinv;
   break;
  case UJ_fomstr:
   alpfun_eval = &UJ_fom_eval;
   wcorinv = &fominv; 
   break;
  case UG_fomstr:
   alpfun_eval = &UG_fom_eval;
   wcorinv = &fominv; 
   break;
  case UQ_fomstr:
   alpfun_eval = &UQ_fom_eval;
   wcorinv = &fominv; 
   break;
  case heidel_unstr:
   alpfun = &Heidel_unstr;
   wcorinv = &unstrinv; 
   break;
  case UJ_equistr:
   alpfun_eval = &UJ_equi_eval;
   wcorinv = &exinv; 
   break;
  case UJ_equimart:
   alpfun_eval = &UJ_equi_eval;
   wcorinv = &exinv; 
   break;
#ifdef EXPERIMENTAL_USERCOR
  case User_cor:
   alpfun = &alpfun_user;
   wcorinv = &wcorinv_user; 
   break;
#endif
  }
 
 //printf("%d = famcode\n",*famcode);
 //printf("%d = corcode\n",*corcode);
 
 double nrec = (double) *n;
 
   matrix Xin = from_S(xin, *n, *p);
   matrix Yin = from_S(yin, *n, 1);
   matrix IDin = from_S(idin, *n, 1);
   matrix Tin = from_S(tin, *n, 1);
   matrix Wall = from_S(weightin, *n, 1);
   matrix OFFin = from_S(OFFSET, *n, 1);
   matrix Binit = from_S(b0, *p, 1);
   matrix Bcur = Binit;
   matrix Alpinit = from_S(alpin, *q, 1);
   double alp = Alpinit.el(0,0);
   int I = cluscount(IDin);
   int cc = Xin.cols();
 //  printf("%d = ncol\n",cc);
 //  printf("%d = I\n",I);
   matrix* X; 
   matrix* Wi; 
   matrix* Y; 
   matrix* TIM; 
   matrix* OFF; 
   X = split(Xin,IDin); /* now X is an I-array of matrices */
   Y = split(Yin,IDin); 
   TIM = split(Tin,IDin); 
   Wi = split(Wall,IDin);
   OFF = split(OFFin,IDin);
 /* starting post-split gee solution */
   double del = 1.0;
   int iter = 0;
   matrix S2 = newmat(*p,*p);
   matrix S2ind = newmat(*p,*p);
   matrix S3 = newmat(*p,1);
   matrix S5 = newmat(*p,*p);
   double phi=0.;
   double u = 0.;
   double duda = 0.;
   double Sum_uut = 0.;
   matrix Mans4;
   int umaxit = 0;
   while (del > *tol & iter < *maxiter ) 
    {
    matrix Mall = apply_elwise(Xin*Bcur+OFFin,(double(*)(double))mu); 
    matrix Rall = Yin - Mall;
    matrix Vall = apply_elwise(Mall,(double(*)(double))variance);
    Vall = mult_like_S(Vall,Wall);
    matrix Aill = apply_elwise(Vall,(double(*)(double))reciproot);
    matrix PRall = mult_like_S(Rall,Aill);
    phi = LZ_scalefun(PRall, IDin, Tin, *p );
 // if fixscale or a UJ_martingale, set phi to 1
    if (*corcode == UJ_equimart | *fixscale == 1) phi = 1.0;
 // if corstr implies an auxiliary estimating function,
 // solve the associated estimating equation
    if (*corcode == UJ_fomstr || *corcode == UG_fomstr || *corcode == UQ_fomstr
              || *corcode == UJ_equistr)
        Mans4 = yua_secant_solve( alpfun_eval, PRall, IDin, Tin, 
                  phi, *p, mat11(alp), *Uagrid_lo, *Uagrid_hi, 
                  *Uagrid_npts, *Uasecant_del, *Uatol, 
                  *Uamaxit);
 // if a martingale est.fun., use martingale resids
    else if (*corcode == UJ_equimart) 
        Mans4 = yua_secant_solve( alpfun_eval, Rall, IDin, Tin, 
                  phi, *p, mat11(alp), *Uagrid_lo, 
                  *Uagrid_hi, *Uagrid_npts, *Uasecant_del, 
                  *Uatol, *Uamaxit);
 // otherwise just evaluate alpfun for this corcode
 // // this is a bad section -- demonstrating that (matrix)q-vector Alpinit
 // is properly handled for fixed alpha
 // TODO: need to handle q>1 and updating of such an alpha
        else Mans4 = (*alpfun)(PRall,IDin, Tin, phi, *p, Alpinit, *Uatol, *Uamaxit );
 
 // presently alpfun can return (alp,u,dudalp,sum.uu')'
    alp = Mans4.el(0,0);
 // patch up the QLS solution for asympt. bias// ChagShults 1999 adj
    if (*corcode == UQ_fomstr) alp = 2.*alp/(1.+alp*alp);  
 
 // for certain corcodes you can get dUda, sum uu' for s.e. alp
    if (*corcode == UJ_fomstr || *corcode == UG_fomstr || *corcode == UJ_equistr)
    {
    u = Mans4.el(1,0);
    duda = Mans4.el(2,0);
    Sum_uut = Mans4.el(3,0);
    //umaxit = (int) Alpout.el(4,0);
    }
 
    iter += 1;
    S2 = 0.*S2;
    S2ind = 0.*S2ind;
    S3 = 0.*S3;
    S5 = 0.*S5;
    qls_val = 0.0;
    ql_val = 0.0;
 if (*verbose) printf("starting iter %d\n",iter);
    for (int i = 0 ; i < I; i++)
      {
      matrix etai = X[i]*Bcur + OFF[i];
      matrix mui = apply_elwise(etai,(double(*)(double))mu);
      matrix dmuideta = apply_elwise(etai,(double(*)(double))dmu_deta);
      matrix Vi = apply_elwise(mui,(double(*)(double))variance);
      Vi = mult_like_S(Vi,Wi[i]);
      matrix oorAi = apply_elwise(Vi,(double(*)(double))reciproot); // 1/sqrt(Ai)
      matrix t1 = mult_like_S(dmuideta,X[i]);
      matrix Di = mult_like_S(oorAi,t1);
      matrix PRi = mult_like_S(oorAi, Y[i]-mui);
 matrix Ri;
 if (*corcode == heidel_unstr) Ri = (matrix)((*wcorinv)(Mans4,Y[i].rows(),TIM[i]));
 // following assumes alpha is scalar!
 // need to clean up design for handing alpha around
      else Ri = (matrix)((*wcorinv)(mat11(alp),Y[i].rows(),TIM[i]));
      S2 = S2 + transp(Di)*Ri*Di;
      S2ind = S2ind + transp(Di)*Di;
      S3 = S3 + transp(Di)*Ri*PRi;
      qls_val = qls_val + (transp(PRi)*Ri*PRi).el(0,0);
      ql_val = ql_val + elsum(dapply_elwise(Y[i],mui,(double(*)(double,double))qlfun));
      S5 = S5 + transp(Di)*Ri*PRi*transp(PRi)*Ri*Di;
      }
    matrix DelB = sweep(S2)*S3;
    del = matmaxabs(DelB);
    Bcur = Bcur + DelB;
    /* if (*verbose) matd(Bcur); */
    }  /* close while del/iter */
  matrix S2i = sweep(S2);
  matrix vs = S2i*S5*S2i;
  to_S(Bcur,bout);
  if (*corcode == heidel_unstr) 
      to_S(Mans4,alpout);
  else to_S(mat11(alp),alpout);
  to_S(mat11(phi),phiout);
  to_S(phi*S2i,bcov_nai);
  to_S(S2i*S5*S2i,bcov_rob);
  *maxiter = iter;
  *Ua = u;
  *dUda = duda;
  *sum_uut = Sum_uut;
 *Uamaxit=umaxit;
 *QLS = qls_val;
 //printf("trace: %lf\n",trace(S5*S2i));
 switch(*famcode)  /* here using p140 of Hardin and Hilbe */
   {
   case gaussian:
     *AIC_PAN = nrec*log(-2.*ql_val/nrec)+2.*trace(S2ind*vs);
     break;
   default:
/*
     printf("-2ql =  %lf \n", -2*ql_val);
     printf("2tr =  %lf \n", 2*trace(S5*S2i));
     printf("2trHH =  %lf \n", 2*trace(S2ind*vs)); */
     *AIC_PAN = -2.*ql_val+2.*trace(S2ind*vs);
     break;
   }
  }
 }
 
 /* Here is end of gee engine */
 
 /*  
 /* extern "C" {
 /* void democc( int* n, double* xin, double* xout )
 /*  {
 /*  matrix dumtim = transp(make_row(0.0,2.0,4.0,6.0));
 /*  matrix Xin = from_S(xin, *n, *n);
 /*  //set_el(Xin,0,0) = 491.3626;
 /*  matrix Z = LDLt(Xin);
 /*  //to_S(form_diag(extract_diag(Z)),xout);
 /*  //to_S(lowutri(Z),xout);
 /*  //matrix Ri = exinv(mat11(.5),3,mat11(0.));
 /*  matrix Ri = fominv(mat11(.5),4,dumtim);
 /*  matrix ZZ = diffL( exinv, mat11(.33333) , 4, dumtim, 0.000001 );
 /* 
 /* GetUjStruc( dB, JJ, BB, fominv, mat11(.33333), 4, dumtim, .000001 )
 /* 
 /*  to_S(JJ,xout);
 /*  }
 /* }  */
 


