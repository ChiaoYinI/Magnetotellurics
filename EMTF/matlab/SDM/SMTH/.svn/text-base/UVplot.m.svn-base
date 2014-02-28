figure
nfmax = 20;
loglog(Sdms.T(1:nfmax),squeeze(Sig_MT(1,1,1:nfmax)),'r')
hold on
loglog(Sdms.T(1:nfmax),squeeze(Sig_MT(2,2,1:nfmax)),'m')
loglog(Sdms.T(1:nfmax),squeeze(Sig_CN(1,1,1:nfmax)),'b')
loglog(Sdms.T(1:nfmax),squeeze(Sig_CN(2,2,1:nfmax)),'c') 
fatlines(gca,2)
legend('MT #1','MT #2','CN #1', 'CN #2')
figure
pol = 1;ri = 'real';
VECplot(T,stcor,pol,ri);

figure
pol = 1;ri = 'real';
VECplot(U,stcor,pol,ri);

