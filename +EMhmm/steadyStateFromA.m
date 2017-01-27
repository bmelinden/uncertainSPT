function [pSS,ok]=steadyStateFromA(A)
numStates=size(A,1);
if(numStates==1)
    pSS=1;
    ok=true;
else
    try
        b=eig(A);
        ok=true;
    catch me
       warning('Steady state could not be computed')
       me
       disp('transition matrix A=')
       disp(nm2str(A))
       pSS=nan(1,size(A,1));
       ok=false;
       return       
    end
    if(sum(abs(b-1)<10*eps)>1)
        warning('Steady state possibly not unique.')
        ok=false;
    end
    bMax=max(b);
    b2nd=max(abs(b(b<bMax)));
    if(abs(bMax-1)<10*eps && ~isempty(b2nd))
        Nss=10*ceil(log(eps)/log(b2nd));
        ASS=(A/bMax)^Nss; % this deals with the possibility that the largest eigenvalue is not exactly 1.
        pSS=ASS(1,:);
    else
        pSS=NaN(1,size(A,1));
        warning('Steady state not found.')
        ok=false;
    end
end
