// note this has to be enriched to allow access to all methods
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
                 friend matrix sweep( matrix );
                 friend matrix ident( int );
                 friend matrix mat11( double );
                 friend double elsum( matrix );

                 void operator=( matrix );  /* assignment */
                 friend void kill( matrix& );

                 };



