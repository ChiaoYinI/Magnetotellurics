function [] = listZ(Z1,Z2,Zhd1,Zhd2,cfile);

fid = fopen(cfile,'w');
nche = Zhd1.nche;
nbt = Zhd1.nbt;
if nche == 3
   comp = { 'Tzx', 'Zxx', 'Zyx'; 'Tzy', 'Zxy', 'Zyy'};
else
   for j = 1:2
      for k = 1:nche
         comp{j,k} = [num2str(k) ',' num2str(j) ];
      end
   end
end
for k = 1:nche
   for j = 1:2
      fprintf(fid,'%s\n\n',['         Component : ',char(comp(j,k))]);

      for ib = 1:nbt
          if Z1.T(ib) > 10
            err1 = real(sqrt(Z1.sigS(j,j,ib)*Z1.sigN(k,k,ib)));
            err2 = real(sqrt(Z2.sigS(j,j,ib)*Z2.sigN(k,k,ib)));
            fprintf(fid,'%8.1f %9.4f %9.4f %10.4f \t %9.4f %9.4f %10.4f\n',...
                Z1.T(ib),real(Z1.z(j,k,ib)),imag(Z1.z(j,k,ib)),err1,...
                real(Z2.z(j,k,ib)),imag(Z2.z(j,k,ib)),err2);
          end
      end
      fprintf(fid,'\n\n','');
   end        
end

