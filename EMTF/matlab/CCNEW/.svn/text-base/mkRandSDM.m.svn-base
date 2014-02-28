function [S] = mkRandSDM(nch,nf,nSamp,V,lambda)
%   makes an array of  nSamp random SDMs each with nch channels
%    and nf samples following the linear model 
%    X_i = sum_k V_k a_ik + e_i
%    Where Cov(a_i) = diag(lambda)
%          Cov(e_i) = I

S = zeros(nch,nch,nSamp)+i*zeros(nch,nch,nSamp);
for j = 1:nSamp
   K = length(lambda);
   A = randn(K,nf)+i*randn(K,nf);
   A = diag(sqrt(lambda/2))*A;
   X = randn(nch,nf)+V*A;
   S(:,:,j) = X*X'/nf;
end
