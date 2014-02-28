function [qcc] = mcSDM(n1,n2,nt,nrep,p,ncc,U,lambda,I1,I2)
% Using monte-carlo calculation comptues quantiles for complex multivariate Gaussian stats
%   qcc: quantiles of canonical coherences
%  ASSUMES: null hypothesis is no coherence
%	USAGE:  [qcc] = mcSDM(n1,n2,nt,nrep,p,ncc,U,lambda,I1,I2);

cc = [];
ncc = min(n1,ncc);ncc = min(n2,ncc);
sqrt2 = sqrt(2);
if(nargin>6)
    N1 = length(I1);N2 = length(I2);
    if(N1 ~= n1 & N2 ~= n2)
        print,'Error: size of I1,I2  do not agree with n1, n2'
        return
    end
    [N,nsig] = size(U);
    Y1 = zeros(n1,nt);
    Y2 = Y1;
    for k=1:nsig
        Y = sqrt(lambda(k))*(randn(1,nt)+i*randn(1,nt));
        Y1 = Y1 + U(I1,k)*Y;
        Y2 = Y2 + U(I2,k)*Y;
    end
else
    nsig = 0;
end
        
for k = 1:nrep
    X1 = randn(n1,nt)+i*randn(n1,nt);
    X2 = randn(n2,nt)+i*randn(n2,nt);
    if nsig > 0
        X1 = X1+Y1;
        X2 = X2+Y2;
    end
    [U1,S1,V1] = svd(X1,0);
    S1 = 1./diag(S1);
    S11hI = U1*diag(S1)*U1';
    [U2,S2,V2] = svd(X2,0);
    S2 = 1./diag(S2);
    S22hI = U2*diag(S2)*U2';
    X12 = S11hI*X1*X2'*S22hI;
    c = svd(X12);
    c = sort(-c);
    c = c(1:ncc).^2;
    cc = [cc c];
end
indP = fix(p*nrep);
qcc = [];
for k = 1:ncc
    C = sort(cc(k,:));
    qcc = [qcc; C(indP)];
end
