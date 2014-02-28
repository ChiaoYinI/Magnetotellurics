function misfit = MsftSmth(nu)

global alpha XX Xd Ruff X d

alpha = Xd/(XX + nu*Ruff);
res = (X*alpha - d);
misfit = sum(sum(conj(res).*res));
