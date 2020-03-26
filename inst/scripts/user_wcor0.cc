#include "MC++.h"
#include "MC++class.h"

matrix alpfun_user(matrix PRin, matrix ID, matrix TIMin, 
       double phi, int p, matrix alpin,
       double atol, int amaxit)
{return mat11(0.0);}

matrix wcorinv_user(matrix alp, int ni, matrix tim) {
// alp, tim included for fixed signature
 return ident(ni);
 }


