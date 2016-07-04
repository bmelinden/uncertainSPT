/*
 * triSym_d1Inv_trjWise.c
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



/* [T0,T1,logDetA]=triSym_triInv_rescale_trjWise(gTot,fTot,trjOne,trjEnd,numTrj)
*
* for each block, defined by g(trjOne(m):trjEnd(m)),
* f(trjOne(m):trjEnd(m)), compute the determinant, and diagonal and first
* off-diagonal of a tridiagonal matrix A with elements
*
* A = [g(1) f(1)    0 0    . . .]
*     [f(1) g(2) f(2) 0    . . .]
*     [   0 f(2) g(3) f(3) . . .]
*     [ . . .              . . .]
*     [ . . .      g(n-1) f(n-1)]
*     [ 0 . . .    f(n-1)  g(n) ]
*
* The matrix is row rescaled by the g-elements, and then transformed to the
* form B = A/diag(g)
* B = [ 1   c(1)  0   0    . . .]
*     [a(2)  1   c(2) 0    . . .]
*     [ 0   a(3)  1  c(3)  . . .]
*     [ . . .              . . .]
*     [ . . .         1   c(n-1)]
*     [ . . .        a(n)    1  ]
* which we then invert and rescale, to get an inverse of the form 
* C = [ T0(1) T1(1)   ?     ?    . . .]
*     [ T1(1) T0(2) T1(2)   ?    . . .]
*     [   ?   T1(2) T0(3) T1(3)  . . .]
*     [ . . .                    . . .]
*     [ . . .         T0(n-1)  T1(n-1)]
*     [ . . .         T1(n-1)   T0(n) ],
*
* and a determinant |A|=|B|*prod(g(j))
*
* Usmani, R. A. (1994). "Inversion of a tridiagonal jacobi matrix".
* Linear Algebra and its Applications. 212-213: 413–414.
* doi:10.1016/0024-3795(94)90414-6
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

/* [T0,T1,logDetA]=triSym_triInv_rescale_trjWise(gTot,fTot,trjOne,trjEnd,numTrj)*/
/* Input Arguments */
#define	G_IN    prhs[0]
#define	F_IN    prhs[1]
#define	TRJONE_IN  prhs[2]
#define	TRJEND_IN  prhs[3]
#define	NUMTRJ_IN  prhs[4]

/* output arguments */
#define	T0_OUT     plhs[0]
#define	T1_OUT      plhs[1]
#define	LOGDETA_OUT     plhs[2]

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]){
    
    /* check number of input/output arguments */
    if (nrhs != 5 )
        mexErrMsgTxt("Five inpit arguments required.");
    if (nlhs != 3)
        mexErrMsgTxt("Three output arguments required.");
    
    /* input parameters */
    double *g,*f,*trjOne,*trjEnd,*numTrj_pt;
    /* output parameters */
    double *T0,*T1,*logDetA;    
    /* temporary variables */
    int gRow,gCol,fRow,fCol,t1Row,t1Col,teRow,teCol; /* matrix dimensions */
    int numTrj,tTot,maxTrjLength,t,nt,tTrj,i;
    double detTrj;
    double *gTrj,*aTrj,*cTrj,*T0trj,*T1trj;    
    mxArray *mx_a,*mx_c,*mx_theta,*mx_phi;
    double *a,*c,*theta,*phi;    
    
    /* retrieve input data */
    /* number of trajectories */
    numTrj_pt=mxGetPr(NUMTRJ_IN);
    numTrj=(int)numTrj_pt[0];    /* inelegant, but works...*/
    
    /* indices to trj starts*/
    t1Row = mxGetM(TRJONE_IN); /* number of rows    */
    t1Col = mxGetN(TRJONE_IN); /* number of columns */
    trjOne=mxGetPr(TRJONE_IN);
    if( t1Row == numTrj && t1Col == 1  ){ /* OK! */
    }else if( t1Row == 1 && t1Col == numTrj ){ /* also OK */
        t1Row=t1Col;
        t1Col=1;
    }else{ /*something weird going on */
        mexErrMsgTxt("trjOne must be numTrj-by-1, or 1-by-numTrj.\n");
    }
    
    /* indices to trj ends*/
    teRow = mxGetM(TRJEND_IN); /* number of rows    */
    teCol = mxGetN(TRJEND_IN); /* number of columns */
    trjEnd=mxGetPr(TRJEND_IN);
    if( teRow == numTrj && teCol == 1  ){ /* OK! */
    }else if( teRow == 1 && teCol == numTrj ){ /* also OK */
        teRow=teCol;
        teCol=1;
    }else{ /*something weird going on */
        mexErrMsgTxt("trjEnd must be numTrj-by-1, or 1-by-numTrj.\n");
    }
    
    /* matrix diagonal */
    gRow = mxGetM(G_IN); /* number of rows    */
    gCol = mxGetN(G_IN); /* number of columns */
    g    =mxGetPr(G_IN);
    if( gRow == trjEnd[numTrj-1] && gCol == 1  ){ /* OK! */
    }else if( gRow == 1 && gCol == trjEnd[numTrj-1] ){ /* also OK */
        gRow=gCol;
        gCol=1;
    }else{ /*something weird going on */
        mexErrMsgTxt("g must be numTrj-by-1, or 1-by-numTrj.\n");
    }
    tTot=gRow; /* total number of data points */
    
    /* matrix off-diagonal */
    fRow = mxGetM(F_IN); /* number of rows    */
    fCol = mxGetN(F_IN); /* number of columns */
    f    =mxGetPr(F_IN);
    if( fRow == tTot && fCol == 1  ){ /* OK! */
    }else if( fRow == 1 && fCol == tTot ){ /* also OK */
        fRow=fCol;
        fCol=1;
    }else{ /*something weird going on */
        mexErrMsgTxt("f must be numTrj-by-1, or 1-by-numTrj.\n");
    }

    /* create output variables */
    T0_OUT = mxCreateDoubleMatrix(gRow, 1, mxREAL);
    T0=mxGetPr(T0_OUT);
    
    T1_OUT = mxCreateDoubleMatrix(gRow, 1, mxREAL);
    T1=mxGetPr(T1_OUT);
    
    LOGDETA_OUT=mxCreateDoubleMatrix(1,1,mxREAL);            
    logDetA=mxGetPr(LOGDETA_OUT);
    
    /* create temporary matrices */
    mx_c=mxCreateDoubleMatrix(gRow, 1, mxREAL);
    mx_a=mxCreateDoubleMatrix(gRow, 1, mxREAL);
    c=mxGetPr(mx_c);
    a=mxGetPr(mx_a);

    /* vectors as long as the longest trajectory */
    maxTrjLength=0;
    for(nt=0;nt<numTrj;nt++){
        tTrj=(int)(trjEnd[nt]-trjOne[nt])+1;
        if(tTrj>maxTrjLength)
            maxTrjLength=tTrj;
    }
    mx_theta=mxCreateDoubleMatrix(maxTrjLength, 1, mxREAL);
    mx_phi  =mxCreateDoubleMatrix(maxTrjLength, 1, mxREAL);
    theta=mxGetPr(mx_theta);
    phi  =mxGetPr(mx_phi);
    
    /* rescale columns: B = A/G; a,c are off-diagonals of B */
    c[tTot-1]=0.0;
    a[0]=0.0;
    logDetA[0]=log(g[tTot-1]);
    for(t=0;t<tTot-1;t++){
        c[t]=f[t]/g[t+1];
        a[t+1]=f[t]/g[t];
        logDetA[0]=logDetA[0]+log(g[t]);
    }
    /* main inversion loop */
    for(nt=0;nt<numTrj;nt++){
        tTrj=(int)(trjEnd[nt]-trjOne[nt])+1; /* length of this trajectory */
        gTrj=&g[(int)trjOne[nt]-1]; /* -1, since trjOne uses matlab index convetion */
        aTrj=&a[(int)trjOne[nt]-1]; /* and these C arrayts do not                   */
        cTrj=&c[(int)trjOne[nt]-1];
        T0trj=&T0[(int)trjOne[nt]-1];
        T1trj=&T1[(int)trjOne[nt]-1];
        
        /* principal minors */
        theta[0]   =1.0;
        theta[1]   =1-aTrj[1]*cTrj[0];
        phi[tTrj-1]=1;
        phi[tTrj-2]=1-aTrj[tTrj-1]*cTrj[tTrj-2];
        for(t=2;t<tTrj;t++){
            theta[t]=theta[t-1]-aTrj[t]*cTrj[t-1]*theta[t-2];
            i=tTrj-1-t;
            phi[i]=phi[i+1]-cTrj[i]*aTrj[i+1]*phi[i+2];
        }
        detTrj=theta[tTrj-1];
        logDetA[0]=logDetA[0]+log(detTrj);
        
        /* diagonal of inverse matrix */
        T0trj[0]=phi[1]/detTrj;
        for(t=1;t<tTrj-1;t++){
            T0trj[t]=theta[t-1]*phi[t+1]/detTrj;
        }
        T0trj[tTrj-1]=theta[tTrj-2]/detTrj;
        
        /* upper off-diagonal of inverse matrix */
        T1trj[0]=-cTrj[0]*phi[2]/detTrj;
        for(t=1;t<tTrj-2;t++){
            T1trj[t]=-cTrj[t]*theta[t-1]*phi[t+2]/detTrj;
        }
        T1trj[tTrj-2]=-cTrj[tTrj-2]*theta[tTrj-3]/detTrj;
        T1trj[tTrj-1]=0.0; /* reset element from earlier trajectories */
        
        
        /* scale back the answer */
        for(t=0;t<tTrj;t++){
            T0trj[t]=T0trj[t]/gTrj[t];
            T1trj[t]=T1trj[t]/gTrj[t];
        }
    }
    /*    for mm=1:length(trjOne)     
    % diagonal
    T0trj(1)=phii(2)/thetaj(tTrj);
    for j=2:tTrj-1
        T0trj(j)=thetaj(j-1)*phii(j+1)/thetaj(tTrj);
    end
    T0trj(tTrj)=thetaj(tTrj-1)/thetaj(tTrj);

    % upper off-diagonal
    T1trj(1)   =-cTrj(1)*phii(3)/thetaj(tTrj);
    for j=2:tTrj-2
        T1trj(j)=-cTrj(j)*thetaj(j-1)*phii(j+2)/thetaj(tTrj);
    end
    T1trj(tTrj-1)=-cTrj(tTrj-1)*thetaj(tTrj-2)/thetaj(tTrj);
    
    % determinant
    logDetA=logDetA+log(thetaj(tTrj));
        
    % scale back rows: inv(A)=inv(G)*inv(B)
    T0trj=T0trj./gTrj;
    T1trj=T1trj./gTrj;

    % reinsert answers
    T0(trjOne(mm):trjEnd(mm))=T0trj;
    T1(trjOne(mm):trjEnd(mm))=T1trj;
end
*/

    /* destroy temporary vriables */
    mxDestroyArray(mx_c);
    mxDestroyArray(mx_a);
    mxDestroyArray(mx_theta);
    mxDestroyArray(mx_phi);
}

