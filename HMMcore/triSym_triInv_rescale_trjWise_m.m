function [T0,T1,logDetA]=triSym_triInv_rescale_trjWise_m(g,f,trjOne,trjEnd,numTrj)
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

% output variables
Ttot=trjEnd(numTrj)-trjOne(1)+1;
T0=zeros(Ttot,1);
T1=zeros(Ttot,1);
logDetA=0; %#ok
% end output variables

% rescale columns: B = A/G;
c=f./g([2:Ttot 1]);
c(Ttot)=0;

a=f./g;
a=a([1 1:Ttot-1]);
a(1)=0;

logDetA=sum(log(g));

for mm=1:length(trjOne)
    gTrj=g(trjOne(mm):trjEnd(mm));
    aTrj=a(trjOne(mm):trjEnd(mm));
    cTrj=c(trjOne(mm):trjEnd(mm));    

    % initialize variables
    nTrj=trjEnd(mm)-trjOne(mm)+1;
    thetaj=zeros(nTrj,1);
    phii  =zeros(nTrj,1);
    T0trj=zeros(nTrj,1);
    T1trj=zeros(nTrj,1);

    % intermediate variables
    thetaj(1)=1;
    thetaj(2)=1-aTrj(2)*cTrj(1);
    phii(nTrj)  =1;
    phii(nTrj-1)=1-aTrj(nTrj)*cTrj(nTrj-1);
    for j=3:nTrj
        thetaj(j)=thetaj(j-1)-aTrj(j)*cTrj(j-1)*thetaj(j-2);
        i=nTrj+1-j;
        phii(i)=phii(i+1)-cTrj(i)*aTrj(i+1)*phii(i+2);
    end    
    % determinant
    logDetA=logDetA+log(thetaj(nTrj));

    % matrix inverse elements
    % diagonal
    T0trj(1)=phii(2)/thetaj(nTrj);
    for j=2:nTrj-1
        T0trj(j)=thetaj(j-1)*phii(j+1)/thetaj(nTrj);
    end
    T0trj(nTrj)=thetaj(nTrj-1)/thetaj(nTrj);

    % upper off-diagonal
    T1trj(1)   =-cTrj(1)*phii(3)/thetaj(nTrj);
    for j=2:nTrj-2
        T1trj(j)=-cTrj(j)*thetaj(j-1)*phii(j+2)/thetaj(nTrj);
    end
    T1trj(nTrj-1)=-cTrj(nTrj-1)*thetaj(nTrj-2)/thetaj(nTrj);
    
    
    % scale back rows: inv(A)=inv(G)*inv(B)
    T0trj=T0trj./gTrj;
    T1trj=T1trj./gTrj;
   
    % reinsert answers
    T0(trjOne(mm):trjEnd(mm))=T0trj;
    T1(trjOne(mm):trjEnd(mm))=T1trj;
end

