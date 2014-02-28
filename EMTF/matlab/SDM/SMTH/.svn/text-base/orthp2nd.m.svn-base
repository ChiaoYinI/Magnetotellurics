function [Q,S,V] = orthp2nd(a,b,N,maxdeg)
%  For an interval [a,b], divided into N+1 equal divisions
%  orthp2nd returns a set of orthogonal polynomials for
%  which the second derivative roughness penalty functional
%  || f'' ||^2 has the simple diagonal form
%  || f'' ||^2 = sum_0^M S(m) c_m^2 where the c_m are coefficients
%   in the expansion of f = sum_0^M c_m Q_m(x).  The polynomials
%   Q_m(x) are returned sampled at the N+1 points in [a,b] in matrix
%   Q;  The columns of V give the coefficients of the polynomials Q_m
 
%N = 100;
%maxdeg = 10;
%a = -2;b = 1;
x = [0:N]/N;
x = a + x*(b-a);
x = x';
P = [];
for k=0:maxdeg
   P = [P x.^k];
end
D = zeros(maxdeg+1,maxdeg+1);
%  D is a matrix which defines derivatives in the 
%   polynomial coefficient domain
for k = 1:maxdeg
   D(k,k+1)=k;
end
%[U,S,V] = svd(P,0);
[Q,R] = qr(P,0);
C = R*D*D/R;
B = C'*C;
[U,S] = eig(B);
S = diag(S);
Q = Q*U;
[S,ind] = sort(S);
Q = Q(:,ind);
V = R\U(:,ind);
