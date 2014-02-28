sysTFfile = 'NIMS_filters.stf';
[SysTF] = rdSysTF(sysTFfile);
figure
amp = abs(1./SysTF.TF);
ph = -(180/pi)*atan2(imag(SysTF.TF),real(SysTF.TF));
freq = SysTF.freq;
loglog(freq,amp(4,:));
set(gca,'FontWeight','demi','FontSize',14)
title('High Pass Filter Response: Amplitude')
figure
semilogx(freq,ph(4,:));
set(gca,'FontWeight','demi','FontSize',14)
title('High Pass Filter Response: Phase')


