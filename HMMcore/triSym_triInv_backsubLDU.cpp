/*
 * triSym_triInv_backsubLDU.c
 * =========================================================================
 *
 * Copyright (C) 2015 Martin Lindén, E-mail: bmelinden@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or any later
 * version.
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * Additional permission under GNU GPL version 3 section 7
 *
 * If you modify this Program, or any covered work, by linking or combining it
 * with Matlab or any Matlab toolbox, the licensors of this Program grant you
 * additional permission to convey the resulting work.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */



/* 
 [y,w0,w1,lnDetT]=triInv_backsub(a,b,x)
 Solve linear system T*y=x, and compute log(det(T)) and the diagonal and
 first off-diagonal entries of inv(T), where T is a symmetric, positive
 definite tridiagonal matrix of the form
 T = [a(1) b(1)    0  . . . . . . ]
     [b(1) a(2) b(2)    0 . . . . ]
     [   0 b(2) a(3) b(3) 0 . . . ]
                  . . . 
    [. . . 0 b(n-2) a(n-1) b(n-1) ]
    [. . . . . . 0  b(n-1)   a(n) ]

T can also be block-diagonal positive definite tridiagonal blocks.

The first two main diagonals of the inverse are computed
inv(T) = [w0(1) w1(1) ...              ]
         [w1(1) w0(2) w1(2) ...        ]
         [...    w1(2) w0(3) w1(3) ... ]
          . . . 
using a slightly modified version of the inversion algorithm from [1],
which should be numerically stable if T is diagonally dominant, i.e., 
a(i)>= -(b(i-1)+b(i)), i=2, ..., n-1, and
a(1)>-b(1), a(n)>-b(n-1).
-------------------------------------------------------
% not computed here (but in the matlab version):
ddMax(1) = max( -[b(i-1)|+b(i)]/a(i)   ), i=2,...,n-1 ,
ddMax(2) = max( -b(n-1)/a_n, b(1)/a(1) ).
-----------------------------------------------------

Finally, an n*1-vector g is returned, which can be used to construct an
LDU-decomposition of T, given by T = L*D^{-1}*L', with  

L = [g(1)    0   . . .            ]
    [b(1) g(2) 0 . . .
                 . . . 
    [           b(n-2) g(n-1)    0]
    [                0 b(n-1) g(n)],

and D=diag(g). This decomposition is used to solve T*y=x by simple
back-substitutions of  (LD^{1})*z=x, and then L'*y=z.

The determinant det(T) = prod(g), i.e., lnDetT = sum(log(g)).
-------------------------------------------------------------------------
1. G. Meurant, “A Review on the Inverse of Symmetric Tridiagonal and Block
Tridiagonal Matrices,” SIAM. J. Matrix Anal. & Appl., vol. 13, no. 3, pp.
707–728, Jul. 1992.  http://dx.doi.org/10.1137/0613045

Modifications: 
1) different index and sign convention for b : 
our b(i) = Meurant's -b(i+1), and Meurant's b(i) = our -b(i-1)
2) We compute the diagonals directly without the detour via the numbers
u,v in [1], to avoid under- and overflow for large systems. From Meurant
[1], but using our definition of b(i), 
one can show that  
w0(1) = 1/d(1);
w0(j) = w0(j-1)*g(j)/d(j);
w1(j) =-w0(j)*b(j)/d(j+1);
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

// [y,w0,w1,lnDetT,ddMax]=triInv_backsub(a,b,x)
// Input Arguments 
#define	A_IN prhs[0]
#define	B_IN prhs[1]
#define	X_IN prhs[2]

// output arguments
#define	Y_OUT      plhs[0]
#define	W0_OUT     plhs[1]
#define	W1_OUT     plhs[2]
#define	LOGDET_OUT plhs[3]

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    
    // check number of input/output arguments 
    if (nrhs != 3 )
        mexErrMsgTxt("Three input arguments required.");
    if (nlhs != 4)
        mexErrMsgTxt("Four output arguments required.");
    
    // input parameters 
    double *a,*b,*x;
    a    =mxGetPr(A_IN);
    b    =mxGetPr(B_IN);
    x    =mxGetPr(X_IN);
    
    // size of input data 
    int aRow,aCol,n;
    aRow = mxGetM(A_IN); // number of rows   
    aCol = mxGetN(A_IN); // number of columns
    n=aRow*aCol;
    
    // output parameters
    double *y,*w0,*w1,*logDet,*ddMax;
    Y_OUT = mxCreateDoubleMatrix(n, 1, mxREAL);
    y=mxGetPr(Y_OUT);
    W0_OUT = mxCreateDoubleMatrix(n, 1, mxREAL);
    w0=mxGetPr(W0_OUT);
    W1_OUT = mxCreateDoubleMatrix(n, 1, mxREAL);
    w1=mxGetPr(W1_OUT);
    LOGDET_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
    logDet=mxGetPr(LOGDET_OUT);
    
    /* compute d, log(|T|), g
      d=zeros(n,1);
      d(n)=a(n);
      for i=n-1:-1:1
           d(i)=a(i)-b(i)^2/d(i+1);
      end
      lnDetT=sum(log(d));
     
      g=zeros(size(d));
      g(1)=a(1);
      for i=2:n
           g(i)=a(i)-b(i-1)^2/g(i-1);
      end */
    
    mxArray *mx_d;
    double *d;
    mx_d=mxCreateDoubleMatrix(n, 1, mxREAL);
    d=mxGetPr(mx_d);
    logDet[0]=0;
    d[n-1]=a[n-1];
    int j;
    logDet[0]=log(d[n-1]);
    for(j=n-2;j>=0;j--){
        d[j]=a[j]-b[j]*b[j]/d[j+1];
        logDet[0]+=log(d[j]);
    }
    mxArray *mx_g;
    double *g;
    mx_g=mxCreateDoubleMatrix(n, 1, mxREAL);    
    g=mxGetPr(mx_g);
    g[0]=a[0];
    for(j=1;j<n;j++){
        g[j]=a[j]-b[j-1]*b[j-1]/g[j-1];
    }
    /*% main diagonals with minimal cancellations:
      w0=zeros(n,1);
      w1=zeros(n,1);
     
      w0(1) =1/d(1);
      for j=2:n
      w0(j)  = w0(j-1)*g(j-1)/d(j);
      w1(j-1)=-w0(j-1)*b(j-1)/d(j);
      end */
    w0[0]=1/d[0];
    for(j=1;j<n;j++){
        w0[j  ]= w0[j-1]*g[j-1]/d[j];
        w1[j-1]=-w0[j-1]*b[j-1]/d[j];
    }
    w1[n-1]=0.0;
    /* % back-substitution
        z=zeros(n,1);
        z(1)=x(1);
        for j=2:n
            z(j)=x(j)-b(j-1)/g(j-1)*z(j-1);
        end
        y=zeros(size(z));
        y(n)=z(n)./g(n);
        for j=n-1:-1:1
            y(j)=(z(j)-b(j)*y(j+1))/g(j);
        end
     */
    mxArray *mx_z;
    double *z;
    mx_z=mxCreateDoubleMatrix(n, 1, mxREAL); 
    z=mxGetPr(mx_z);
    z[0]=x[0];
    for(j=1;j<n;j++){
        z[j]=x[j]-b[j-1]/g[j-1]*z[j-1];
    }
    y[n-1]=z[n-1]/g[n-1];
    for(j=n-2;j>=0;j--){
        y[j]=(z[j]-b[j]*y[j+1])/g[j];
    }
    
    /* % check diagonal stability: skip!
        ddMax(1)=max(-(b(2:n)+b(1:n-1))./a(2:n));
        ddMax(2)=max(-b(1)/a(1),-b(n-1)/a(n));
     */
    

    /* destroy temporary vriables */
    mxDestroyArray(mx_d);
    mxDestroyArray(mx_g);
    mxDestroyArray(mx_z);
}
