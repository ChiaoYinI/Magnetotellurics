function [Z,Zhd] = Zstruc(cfile);
% Usage : [Z,Zhd] = Zstruc(cfile);

[z,sig_s,sig_e,periods,ndf,stdec,orient,nch,nche,nbt,chid,csta] = Z_in(cfile);

Z = struct('T',periods,'nf',ndf,'z',z,'sigS',sig_s,'sigN',sig_e);
Zhd = struct('nbt',nbt,'nch',nch,'nche',nche,...
	'decl',stdec,'orient',orient,'chid',chid,'csta',csta);
return