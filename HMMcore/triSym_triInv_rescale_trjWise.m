function [T0,T1,logDetA]=triSym_triInv_rescale_trjWise(g,f,trjOne,trjEnd,numTrj)
% [T0,T1,logDetA]=triSym_triInv_rescale_trjWise(gTot,fTot,trjOne,trjEnd,numTrj)
%
% for each block, defined by g(trjOne(m):trjEnd(m)),
% f(trjOne(m):trjEnd(m)), compute the determinant, and diagonal and first
% off-diagonal of a tridiagonal matrix A with elements. 
%
% A = [g(1) f(1)    0 0    . . .]
%     [f(1) g(2) f(2) 0    . . .]
%     [   0 f(2) g(3) f(3) . . .]
%     [ . . .              . . .]
%     [ . . .      g(n-1) f(n-1)]
%     [ 0 . . .    f(n-1)  g(n) ]
%
% The matrix is row rescaled by the g-elements, and then transformed to the
% form B = A/diag(g)
% B = [ 1   c(1)  0   0    . . .]
%     [a(2)  1   c(2) 0    . . .]
%     [ 0   a(3)  1  c(3)  . . .]
%     [ . . .              . . .]
%     [ . . .         1   c(n-1)]
%     [ . . .        a(n)    1  ]
% which we then invert and rescale, to get an inverse of the form 
% C = [ T0(1) T1(1)   ?     ?    . . .]
%     [ T1(1) T0(2) T1(2)   ?    . . .]
%     [   ?   T1(2) T0(3) T1(3)  . . .]
%     [ . . .                    . . .]
%     [ . . .         T0(n-1)  T1(n-1)]
%     [ . . .         T1(n-1)   T0(n) ],
%
% and a determinant |A|=|B|*prod(g(j))
%
% Usmani, R. A. (1994). "Inversion of a tridiagonal jacobi matrix".
% Linear Algebra and its Applications. 212-213: 413–414.
% doi:10.1016/0024-3795(94)90414-6
%
% ML 2014-11-10
