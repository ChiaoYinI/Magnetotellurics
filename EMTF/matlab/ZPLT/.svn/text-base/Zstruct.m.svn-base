function [Z,Zhd] = Zstruct(cfile);
% Usage : [Z,Zhd] = Zstruct(cfile);

[z,sig_s,sig_e,periods,ndf,stdec,orient,nch,nche,nbt,chid,csta] = Z_in(cfile);
z = reshape(z,[2,nche,nbt]);
sig_s = reshape(sig_s,[2,2,nbt]);
sig_e = reshape(sig_e,[nche,nche,nbt]);

Z = struct('T',periods,'nf',ndf,'z',z,'sigS',sig_s,'sigN',sig_e);
Zhd = struct('nbt',nbt,'nch',nch,'nche',nche,...
	'decl',stdec,'orient',orient,'chid',chid,'csta',csta);
return
