function [V,A,T,B1,B2] = CC_TF(S,I1,i2,B1,B2,JOB)
%  Usage: [V,A,T,B1,B2] = CC_TF(S,I1,I2,B1,B2,JOB);

B1half = sqrtm(B1);
B2half = sqrtm(B2);
S21 = S(I2,I1);
S21t = B2half\S(I2,I1)/B1half';
S11t = B1half\S(I1,I1)/B1half;
S22t = B2half\S(I2,I2)/B2half;
[U1,s,U2] = svd(S21t,0);

nCoh = sum(s > JOB.sig1);
U1 = U1(:,1:nCoh);
U2 = U2(:,1:nCoh);
Shat = [ U1'*S11t*U1 U1'*S21t'*U2; ...
         U2'*S21t*U1 U2'*S22t*U2];
[W,LAMBDA] = eig(Shat);
[temp,ind] = sort(diag(LAMBDA));
nt = length(ind);
ind = ind([nt:-1:1]);
W = W(:,ind(1:nCoh);
LAMBDA = temp([nt:-1:nt-nCoh+1]) ;

T = W(I1,:)/W(I2,:);
A = W(I1,:)*(diag(LAMBDA-1))*(W(I1,:)');
U1 = B1half*U1;
U2 = B2half*U2;

%Shat = [U1*A*U1'      U1*A*U2'   ; ...
%        U2*T*A*U1   U2*T*A*T*U2'];
B1 = S(I1,I1)-U1*A*U1';
B2 = S(I2,I2)-U2*T*A*T*U2'];
