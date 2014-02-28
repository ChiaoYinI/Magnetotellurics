function [Zall,Z2x2all] = ZallRead(fileList)

%    Usage:  [Zall,Z2x2all] = ZallRead(fileList);
%   reads in files in File list, makes cell arrays of impedance
%    data structures
%   NOTE:  Zall is the "raw" TF structure, containing all TFs
%    between two input channels (generally Hx, Hy) and all other channels.
%    With multiple dipoles (EMAP, or parallel dipoles) all are included
%    in this structure, with TFs expressed in the raw measurement coordinates
%    Z2x2all "pairs up" Ex,Ey to form one or more impedandes, and converts
%    these into a standard orthogonal coordinate system.

%   need to use the correct version of Z_in!
%    path('/home/lenz/data/PKDSAO/matlab/EMTF/IN',path);
    nFiles = length(fileList);
    k = 0;
    NBT = zeros(nFiles,1);
    for iFile = 1:nFiles
        cfile = fileList{iFile};
       k = k +1;

       % see if the file exists 
      if ( fopen(cfile) > 0 )
         %fclose(fopen(cfile));
         %18Aug05, Too many files error?  change above to below
         fclose('all');
         %if so, call Z_in to get the impedances from the Z-file
	 [z,sig_s,sig_e,T,ndf,stdec,orient, ...
		nch,nche,nbt,chid,csta] = Z_in(cfile);%display 'Z_in ok';cfile
         %  put into temporary data structure
         NBT(k) = nbt;
         if (nbt > 0) 
           Ztemp = struct('type','Zraw','miss',0,'good',1, ...
             'Z',z,'S',sig_s,'E',sig_e,'T',T,'ndf',ndf,...
             'orient',orient,'nch',nch,'nche',nche,'nbt',nbt,...
             'chid',chid,'csta',csta);

           %check that ecomp will run properly... does E field data exist?
           %fubar = ecompworkskk(nche,chid)
           %if fubar ==0
               [xy,yx,hz,xypairs,DipoleSetup]=ecomp(nche,chid);
           %end;
           %   convert to impedances in geographic coordinates
           
           [Z2x2,SIG_S,SIG_E] = z_to_imp(z,sig_e,sig_s,nche,xypairs,orient);
           [dum,nImpxNbt] = size(Z2x2);
           nImp = nImpxNbt/nbt;
           temp = reshape(Z2x2,[2,2,nbt,nImp]);
           for ll = 1,nbt;
              for kk = 1,nImp;
                  temp(:,:,ll,kk) = squeeze(temp(:,:,ll,kk)).';
              end
           end
           Z2x2 = temp;
           temp = reshape(SIG_S,[2,2,nbt,nImp]);
           %   should this be transposed also????
           SIG_S = temp;
           temp = reshape(SIG_E,[2,2,nbt,nImp]);
           for ll = 1,nbt;
              for kk = 1,nImp;
                temp(:,:,ll,kk) = squeeze(temp(:,:,ll,kk)).';
              end
           end
           SIG_E = temp;

           %   and add to cell array of time segment impedance structures
           Z2x2temp = struct('type','Z2x2','miss',0,'good',1,...
                 'Z',Z2x2,'S',SIG_S,'E',SIG_E,'T',T,'ixy',xypairs);
           Zall{k} = Ztemp;
           Z2x2all{k} = Z2x2temp;
	   clear Ztemp ;
	   clear Z2x2temp ;
         else    %  nbt = 0
           Zall{k} = struct('miss',1,'good',0);
           Z2x2all{k} = struct('miss',1,'good',0);
         end    %  end of nbt > 0
      else    %  file not found ... mark as missing
         Zall{k} = struct('miss',1,'good',0);
         Z2x2all{k} = struct('miss',1,'good',0);
      end		%   End of if file exists
   end		%  End of loop over files
   %  now check to see that all have same number of bands
   NBTnormal = median(NBT(NBT~=0));
   for k = 1:nFiles
      if(NBT(k) ~= NBTnormal)
         Z2x2all{k}.miss = 1;
         Z2x2all{k}.good = 0;
      end
   end
