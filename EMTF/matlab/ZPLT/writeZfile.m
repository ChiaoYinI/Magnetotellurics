function [] = writeZfile(Z,cfile)
 
% Usage:    [] = writeZfile(Z,cfile);
% Writes impedance structure Z (as obtained from readZfile) to file cfile

rnan = -999.;
cform100 = 'station    :%20s\n';
cform105 = 'coordinate %9.3f %9.3f  declination %8.2f\n';
cform110 = 'number of channels %3d   number of frequencies %4d\n';
cform115 = '%5d %8.2f %8.2f %3s %6s\n';
cform120 = ['period : %12.5f    decimation level %3d'...
        	'    freq. band from %4d to %4d\n'];
cform125 = 'number of data point %6d sampling freq. %7.3f Hz\n';
cform140 = ['%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e'...
           '%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e%12.4e'];
 
fid =  fopen(cfile,'w');
temp = ' **** IMPEDANCE IN MEASUREMENT COORDINATES ****';
fprintf(fid,'%s\n',temp);
temp = ' ********** WITH FULL ERROR COVARINCE**********';
fprintf(fid,'%s\n',temp);
fprintf(fid,'%s\n',Z.chead);
fprintf(fid,cform100,Z.stname);
fprintf(fid,cform105,Z.stcor,Z.theta0);
fprintf(fid,cform110,Z.Nch,Z.nbt);
fprintf(fid,' orientations and tilts of each channel \n');
for k = 1:Z.Nch
   fprintf(fid,cform115,k,Z.orient(:,k),char(Z.sta{k}),char(Z.chid{k}));
end
fprintf(fid,'\n')

for ib = 1:Z.nbt
   fprintf(fid,cform120,Z.T(ib),Z.level(ib),Z.ibandlim(:,ib));
   fprintf(fid,cform125,Z.ndf(ib),Z.sampRate(ib));
   fprintf(fid,' Transfer Functions\n');
   for k = 1:Z.Nche
       temp = Z.TF(:,k,ib);
       temp = [real(temp) imag(temp)]';
       temp(isnan(temp)) = rnan;
       fprintf(fid,cform140,temp);
       fprintf(fid,'\n');
   end
   fprintf(fid,' Inverse Coherent Signal Power Matrix\n');
   for k = 1:2
       temp = Z.SIG_S(1:k,k,ib);
       temp = [real(temp) imag(temp)]';
       temp(isnan(temp)) = rnan;
       fprintf(fid,cform140,temp);
       fprintf(fid,'\n');
   end
   fprintf(fid,' Residual Covariance\n');
   for k = 1:Z.Nche
       temp = Z.SIG_E(1:k,k,ib);
       temp = [real(temp) imag(temp)]';
       temp(isnan(temp)) = rnan;
       fprintf(fid,cform140,temp);
       fprintf(fid,'\n');
   end
end 
fclose(fid);
