function y=backsubLDU(g,b,x)
% y=backsubLDU(g,b,x)
% solve tridoagonal system on the form T*y=x, where the symmetric
% tridiagonal matrix T has been decomposed to and LDU form T = L*D^{-1}*L',
% with 
%
% L = [g(1)    0   . . .            ]
%     [b(1) g(2) 0 . . .
%                  . . . 
%     [           b(n-2) g(n-1)    0]
%     [                0 b(n-1) g(n)],
%
% and D=diag(g).
% y is computed by 2step backsubstitution, by solving first (LD^{1})*z=x,
% and then L'*y=z.
T=length(g);
z=zeros(T,1);
z(1)=x(1);
for j=2:T
    z(j)=x(j)-b(j-1)/g(j-1)*z(j-1);
end
y=zeros(size(z));
y(T)=z(T)./g(T);
for j=T-1:-1:1
    y(j)=(z(j)-b(j)*y(j+1))/g(j);
end
    