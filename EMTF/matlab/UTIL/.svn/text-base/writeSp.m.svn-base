function [status] = writeSp(metaData,fileRoot)

%  Usage:  writeSp(metaData,fileRoot)
%
%  writeSp should use metaData.SystemResponse, if present;
%  an example for NIMS may be created with SS/NIMsysRsp.m


%   First some "data" to write to all  of the files
ChID = {'Hx','Hy','Hz','Ex','Ey'};

%   NOTE we use the same filters for all sites; slightly different time delays
%    for each magnetic channel.  These numbers were correct for 8 hz; might not be
%    exactly correct for 1 hz
MagFilters = {'0.01  2  conversion nT/count, # of filters',...
   'PZ  :3 poles, Butterworth low-pass filter',...
   '0 3 1984.31   (-6.28319,10.8825) (-6.28319,-10.8825) (-12.5664,0.0)',...
   'DT time offset (seconds)',...
   '.2455','.2365','.2275'};
ElecFilters = {'2.44141221047903e-06  3  conversion mV/count, # of filters', ...
   'PZ  :1 pole, Butterworth high-pass filter', ...
   '1 1 1  (0.0,0.0)  (-1.666670E-04,0.0)',...
   'PZ  :5 poles, Butterworth low-pass filter',...
   '0 5 313384   (-3.88301,11.9519) (-3.88301,-11.9519) (-10.1662,7.38651) (-10.1662,-7.38651) (-12.5664,0.0)',...
   'DT time offset (seconds)',...
   '.1525'};

lat = metaData.lat;
lon = metaData.lon;
decl = metaData.decl;
elev = metaData.elev;
gain = metaData.gain_char;
ExAz = metaData.Ex_wire_azimuth;
EyAz = metaData.Ey_wire_azimuth;
ExLength = metaData.Ex_wire_length;
EyLength = metaData.Ey_wire_length;
MagAz = [decl decl+90. 0. ];
ElecAz = [ExAz EyAz]+decl;
DipoleLengths = [ExLength EyLength]/1000;

if nargin < 2
    fileRoot = metaData.runID;
end

%  open output sp file
fileName = [fileRoot '.sp'];
fid = fopen(fileName,'w');

%  site ID
fprintf(fid,'%s\n',metaData.runID);

%  lat, lon, elev
fprintf(fid,'%10.4f %10.4f %7.0f\n',[lat,lon,elev]);

%  coordinate system orientation (express everything in geographic)
fprintf(fid,'%6.2f\n',0.0);

%  nch  ... assume always 5!
fprintf(fid,'%d\n',5);

%  sampling rate: delta t in seconds (need to input, or put into metadata)
fprintf(fid,'%10.4f\n',metaData.dt);

%  clock offset, linear drift (fossil code!)
fprintf(fid,'%d %d \n',[0,0]);

%  Now loop over channels:
%    First 3 mag channels
for ich = 1:3
   fprintf(fid,'%s \n',ChID{ich});
   fprintf(fid,'%8.2f %8.2f\n',[MagAz(ich) 0.0]);
   for k = 1:4
      fprintf(fid,'%s\n',MagFilters{k});
   end
   fprintf(fid,'%s\n',MagFilters{4+ich});
end

% Gain always H(igh): software gain already applied where necessary
value(1:2) = 1;

% Then 2 electric channels
for ich = 1:2
   fprintf(fid,'%s \n',ChID{ich+3});
   fprintf(fid,'%8.2f %8.2f %d %f\n',[DipoleLengths(ich),ElecAz(ich),0,value(ich)]);
   for k = 1:7 
      fprintf(fid,'%s\n',ElecFilters{k});
   end
end   

%  close file
status = fclose(fid);
