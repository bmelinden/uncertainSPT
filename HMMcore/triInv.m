function [w0,w1,ddMax,g,lnDetT]=triInv(a,b)
% [w0,w1,ddMax,g,lnDetT]=triInv(a,b)
% computes first two principal diagonals of the inverse of a symmetric,
% positive definite tridiagonal matrix defined through
% T = [a(1) b(1)    0  . . . . . . ]
%     [b(1) a(2) b(2)    0 . . . . ]
%     [   0 b(2) a(3) b(3) 0 . . . ]
%                  . . . 
%     [. . . 0 b(n-2) a(n-1) b(n-1) ]
%     [. . . . . . 0  b(n-1)   a(n) ]
%
% T can also be block-diagonal positive definite tridiagonal blocks.
%
% The first two main diagonals of the inverse are computed
% inv(T) = [w0(1) w1(1) ...              ]
%          [w1(1) w0(2) w1(2) ...        ]
%          [...    w1(2) w0(3) w1(3) ... ]
%           . . . 
% using a slightly modified version of the inversion algorithm from [1],
% which should be numerically stable if T is diagonally dominant, i.e., 
% a(i)>= -(b(i-1)+b(i)), i=2, ..., n-1, and
% a(1)>-b(1), a(n)>-b(n-1).
%
% ddMax(1) = max( -[b(i-1)|+b(i)]/a(i)   ), i=2,...,n-1 ,
% ddMax(2) = max( -b(n-1)/a_n, b(1)/a(1) ).
%
% Finally, an n*1-vector g is returned, which can be used to construct an
% LDU-decomposition of T, given by T = L*D^{-1}*L', with  
%
% L = [g(1)    0   . . .            ]
%     [b(1) g(2) 0 . . .
%                  . . . 
%     [           b(n-2) g(n-1)    0]
%     [                0 b(n-1) g(n)],
%
% and D=diag(g).
% 
% The determinant det(T) = prod(g), i.e., lnDetT = sum(log(g)).
% 
% 1. G. Meurant, “A Review on the Inverse of Symmetric Tridiagonal and Block
% Tridiagonal Matrices,” SIAM. J. Matrix Anal. & Appl., vol. 13, no. 3, pp.
% 707–728, Jul. 1992.  http://dx.doi.org/10.1137/0613045
%
% Modifications: 
% 1) different index and sign convention for b : 
% our b(i) = Meurant's -b(i+1), and Meurant's b(i) = our -b(i-1)
%
% 2) We compute the diagonals directly without the detour via the numbers
% u,v in [1], to avoid under- and overflow for large systems. From Meurant
% [1], but using our definition of b(i), 
% one can show that  
% w0(1) = 1/d(1);
% w0(j) = w0(j-1)*g(j)/d(j);
% w1(j) =-w0(j)*b(j)/d(j+1);

% check diagonal stability
n=length(a);
ddMax(1)=max(-(b(2:n)+b(1:n-1))./a(2:n));
ddMax(2)=max(-b(1)/a(1),-b(n-1)/a(n));

% compute d
d=zeros(n,1);
d(n)=a(n);
for i=n-1:-1:1
    d(i)=a(i)-b(i)^2/d(i+1);
end

% ln(|T|)
lnDetT=sum(log(d));

% compute g(i)=delta(i)
g=zeros(size(d));
g(1)=a(1);
for i=2:n
    g(i)=a(i)-b(i-1)^2/g(i-1);
end
% check: OK!
%check_lnDet=(lnDetT-sum(log(g)))/lnDetT

% main diagonals with minimal cancellations
w0=zeros(n,1);
w1=zeros(n,1);

w0(1) =1/d(1);
for j=2:n
    w0(j)  = w0(j-1)*g(j-1)/d(j);
    w1(j-1)=-w0(j-1)*b(j-1)/d(j);
end



